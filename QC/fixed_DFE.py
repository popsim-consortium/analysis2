"""
Simulate a single region with a genomic element of a given length. Flanking
regions can be included.
"""

import numpy as np
import stdpopsim
import argparse
import sys

import moments  # can be installed via `conda install moments`

import matplotlib.pylab as plt

import warnings

# I know this probably isn't good, but my god so many annoying warnings
warnings.filterwarnings("ignore")


def make_parser():
    ADHF = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser("flank_simulation.py", formatter_class=ADHF)
    parser.add_argument("--seed", required=True, type=int, help="Random seed.")
    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        "--replicates",
        type=int,
        default=1,
        help="Number of replicates to run, taking sum of resulting spectra.",
    )
    optional.add_argument(
        "--population_size",
        "-N",
        type=int,
        default=10000,
        help="Diploid population size, defaults to 10,000.",
    )
    optional.add_argument(
        "--selection_coefficient",
        "-s",
        type=float,
        default=0.0,
        help="Selection coefficient for selected mutations, defaults to 0.",
    )
    optional.add_argument(
        "--dominance_coefficient",
        type=float,
        default=0.5,
        help="Dominance coefficient for selected mutations, defaults to 0.5.",
    )
    optional.add_argument(
        "--fraction_selected",
        "-f",
        type=float,
        default=0.7,
        help="Fraction of new mutations in selected regions that are selected.",
    )
    optional.add_argument(
        "--num_samples",
        "-n",
        type=int,
        default=100,
        help="Number of (haploid) samples.",
    )
    optional.add_argument(
        "--element_length",
        "-L",
        type=int,
        default=100000,
        help="Length of the genomic element.",
    )
    optional.add_argument(
        "--flanking_length",
        type=int,
        default=0,
        help="Length of the flanking regions.",
    )
    optional.add_argument(
        "--burn_in", type=int, default=10, help="Burn in factor to pass to slim.",
    )
    optional.add_argument(
        "--scaling_factor",
        type=int,
        default=10,
        help="Scaling factor to pass to slim.",
    )
    optional.add_argument(
        "--plotting",
        type=str,
        default=None,
        help="Plotting output file name. If none given, do not plot.",
    )
    return parser


def DFE(selCoeff, fracSelected, domCoeff=0.5):
    """
    Following the DFE construction from the README in the analysis2 repo.
    
    selCoeff: the (unscaled) selection coefficient
    fracSelected: the proportion of new mutations in selected element that are
        selected. 1 - fracSelected will be neutral.
    """
    neutral = stdpopsim.ext.MutationType()
    h = domCoeff
    selected = stdpopsim.ext.MutationType(
        dominance_coeff=domCoeff, distribution_type="f", distribution_args=[selCoeff]
    )
    return {
        "mutation_types": [neutral, selected],
        "proportions": [1 - fracSelected, fracSelected],
    }


def set_up_contig(args, species):
    """
    Sets up the contig with the selected element and flanking region.
    """
    F = args.flanking_length
    L = args.element_length
    total_length = 2 * F + L
    contig = species.get_contig(length=total_length)
    # need to clear genome mutation types before adding them ???
    contig.clear_genomic_mutation_types()
    intervals = np.array([[F, F + L]])
    contig.add_genomic_element_type(
        intervals=intervals,
        **DFE(
            args.selection_coefficient,
            args.fraction_selected,
            domCoeff=args.dominance_coefficient,
        )
    )
    return contig


def get_spectra(ts):
    sel = np.zeros(ts.num_samples + 1)
    neu = np.zeros(ts.num_samples + 1)
    for tree in ts.trees():
        for mut in tree.mutations():
            if len(mut.metadata["mutation_list"]) != 1:
                continue
            s = mut.metadata["mutation_list"][0]["selection_coeff"]
            i = tree.get_num_samples(mut.node)
            if s == 0:
                neu[i] += 1
            else:
                sel[i] += 1
    return sel, neu


def get_expected_spectra(args, contig):
    """
    Get the expected neutral and selected frequency spectra, given simulation parameters
    and assuming sites are unlinked.

    In moments and dadi, scaled selection coefficients are defined as S=2*Ne*s, since if
    A is the derived allele and a is the ancestral allele with aa genotypes having
    relative fitness 1, then Aa have fitness 1+s and AA have fitness 1+2s.

    Slim appears to set the fitess of Aa as 1+s/2 and AA as 1+s. This needs to be double
    checked, but here I've set S=Ne*s instead of 2*Ne*s to get a better match between
    slim simulations and expectations using moments.
    """
    theta = 4 * args.population_size * contig.mutation_rate * args.element_length
    S = args.population_size * args.selection_coefficient
    expected_sel = (
        moments.LinearSystem_1D.steady_state_1D(
            args.num_samples, gamma=S, h=args.dominance_coefficient
        )
        * theta
        * args.replicates
        * args.fraction_selected
    )
    expected_neu = (
        moments.LinearSystem_1D.steady_state_1D(args.num_samples)
        * theta
        * args.replicates
        * (1 - args.fraction_selected)
    )
    return moments.Spectrum(expected_sel), moments.Spectrum(expected_neu)


def plot_spectra_comparison(args, simulation, expected):
    fs_sel, fs_neu = simulation
    exp_sel, exp_neu = expected
    fig = plt.figure(1, figsize=(8, 4))
    fig.clf()
    ax1 = plt.subplot(1, 2, 1)
    ax1.semilogy(fs_sel, "-o", ms=5, lw=1, mfc="w", label="Slim")
    ax1.semilogy(exp_sel, "-o", ms=2, lw=1, label="Moments")
    ax1.set_ylabel("Count")
    ax1.set_xlabel("Allele frequency")
    ax1.legend()
    ax1.set_title("Selected mutations")
    ax2 = plt.subplot(1, 2, 2)
    ax2.semilogy(fs_neu, "-o", ms=6, lw=1, mfc="w", label="Slim")
    ax2.semilogy(exp_neu, "-o", ms=3, lw=1, label="Moments")
    ax2.set_xlabel("Allele frequency")
    ax2.set_title("Nuetral mutations")
    fig.tight_layout()
    plt.savefig(args.plotting)

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    species = stdpopsim.get_species("HomSap")
    model = stdpopsim.PiecewiseConstantSize(args.population_size)
    samples = model.get_samples(args.num_samples)
    contig = set_up_contig(args, species)

    engine = stdpopsim.get_engine("slim")

    np.random.seed(args.seed)
    seeds = np.random.randint(0, np.iinfo(np.uint32).max, args.replicates)

    fs_sel = np.zeros(args.num_samples + 1)
    fs_neu = np.zeros(args.num_samples + 1)

    for ii, random_seed in enumerate(seeds):
        print("running replicate", ii + 1, "of", args.replicates)
        ts = engine.simulate(
            model,
            contig,
            samples,
            seed=random_seed,
            slim_scaling_factor=args.scaling_factor,
            slim_burn_in=args.burn_in,
        )

        sel, neu = get_spectra(ts)
        fs_sel += sel
        fs_neu += neu

    fs_sel = moments.Spectrum(fs_sel)
    fs_neu = moments.Spectrum(fs_neu)
    exp_sel, exp_neu = get_expected_spectra(args, contig)

    if args.plotting is not None:
        plot_spectra_comparison(args, (fs_sel, fs_neu), (exp_sel, exp_neu))
