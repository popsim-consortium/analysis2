"""
Code for generating plots.
"""
import pandas as pd
from math import sqrt
from pathlib import Path
import glob
import os
import sys
import matplotlib
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import subprocess

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
#sns.set_style("darkgrid")
COLOURS = sns.color_palette("colorblind")

def plot_sfs(s, outfile):
    """
    Plot the SFS for this simulation
    """
    bins = [n + 1 for n in range(len(s[0]))]
    vals = []
    for i in range(len(s)):
        vals.append([int(x) for x in s[i]])
    if len(s) == 2:
        f, ax = plt.subplots(
            1, 2, sharey=True, tight_layout=True, figsize=(8, 3))
        ax[0].bar(bins, vals[0])
        ax[1].bar(bins, vals[1])
        ax[0].set_title("all sites")
        ax[1].set_title("masked sites")
        ax[0].set_ylabel("counts")
        ax[0].set_xlabel("derived allele frequency")
        ax[1].set_xlabel("derived allele frequency")
        f.savefig(outfile, bbox_inches='tight')
        plt.close()
    else:
        f, ax = plt.subplots(figsize=(3, 3))
        ax.bar(bins, vals[0])
        ax.set_title("all sites")
        ax.set_ylabel("counts")
        ax.set_xlabel("derived allele frequency")
        f.savefig(outfile, bbox_inches='tight')
        plt.close()


def gather_inference_results(output_dir, demog, output, method, chrm_mask, 
                             annot_mask, pops_size_dict, slim_scaling_factor, gen_time):
    infiles = glob.glob(f"{output_dir}/inference/{demog}/{method}/**/{method}_estimated_Ne.*", recursive=True)
    header = "method,population,nsamp,DFE,annotations,year,Ne,seed,chrm_mask,annot_mask,slim_scaling_factor"
    if method == "stairwayplot":
        header += ",Ne_02_5,Ne_97_5"
    elif method == "msmc":
        header += ",n_genomes"
    elif method == "gone":
        header += ",G_value"
    elif method == "smcpp":
        header += ",n_genomes"
    with open(output, 'w') as f:
        f.write(f"{header}\n")
        for infile in infiles:
            infile_path = Path(infile)
            infile_parts = infile_path.parts
            pop = infile_parts[-2]
            size = pops_size_dict[pop]
            seed = infile_parts[-3]
            annot = infile_parts[-4]
            dfe = infile_parts[-5]
            if chrm_mask is None:
                chrm_mask_i = "none"
            else:
                chrm_mask_i = Path(chrm_mask).name
            if annot_mask == 'none':
                annot_mask_i = 'none'
            else:
                annot_mask_i = annot
            if method == "stairwayplot":
                nt = pd.read_csv(infile, sep="\t", skiprows=5)
                nt.columns = nt.columns.str.replace('[%,.]','')
                for row in nt.itertuples():
                    f.write(f'{method},{pop},{size},{dfe},{annot},{row.year},{row.Ne_median},{seed},{chrm_mask_i},{annot_mask_i},{slim_scaling_factor},{getattr(row, "Ne_25")},{getattr(row, "Ne_975")}\n')
            elif method == "msmc":
                nt = pd.read_csv(infile, sep="\t")
                for row in nt.itertuples():
                    f.write(f'{method},{pop},{size},{dfe},{annot},{row.year},{row.Ne},{seed},{chrm_mask_i},{annot_mask_i},{slim_scaling_factor},{row.n_samp}\n')
            elif method == "smcpp":
                nt = pd.read_csv(infile, sep=",")
                for row in nt.itertuples():  #row[1], getattr(row, "name"), row.name
                    f.write(f'{method},{pop},{size},{dfe},{annot},{row.x},{row.y},{seed},{chrm_mask_i},{annot_mask_i},{slim_scaling_factor},{2}\n')
            elif method == "gone":
                q_file = infile_path.parent / "OUTPUT_gone"
                with open(q_file) as q:
                    i=1
                    for lin in q:
                        if "Phase" in lin:
                            break
                        i += 1
                ld_pairs = pd.read_csv(q_file, sep=" ", skiprows=i+2, names=["ld_pairs","avg_c", "avg_d2", "gens"])
                nt = pd.read_csv(infile, sep="\t", skiprows=1)
                nt = nt[nt["Generation"] <= 200]
                for row in nt.itertuples():
                    generation = row.Generation
                    Ne = row.Geometric_mean
                    ld_bin = ld_pairs[ld_pairs["gens"] <= generation]
                    while len(ld_bin.index) == 0:
                        generation += 1
                        ld_bin = ld_pairs[ld_pairs["gens"] <= generation]
                    G_val = (size * sqrt(ld_bin.iloc[-1]["ld_pairs"])) / Ne
                    f.write(f'{method},{pop},{size},{dfe},{annot},{row.Generation * gen_time},{Ne},{seed},{chrm_mask_i},{annot_mask_i},{slim_scaling_factor},{G_val}\n')
            else:
                print("Error: Method not recognized")
                sys.exit(1)


def popn_coal_rate(model, pop_id, pops_dict, generation_time, steps):
    """
    returns tuple (coal_rate, P, steps) for pop_id
    conditional on the model and configuration

    coal_rate, P, steps = popn_coal_rate(
        model, pop_id, n_samp, generation_time)
    """
    ddb = model.model.debug()
    if steps is None:
        if(len(ddb.epochs) == 1):
            end_time = 100000
        else:
            end_time = int(ddb.epochs[-2].end_time) + 10000
        steps = np.linspace(1, end_time, 1000)
    sample_dict = {}
    for pop in pops_dict:
        size = 0
        if pop == pop_id:
            size=pops_dict[pop]
        sample_dict[pop] = size
    coal_rate, P = ddb.coalescence_rate_trajectory(steps=steps,
                                                   lineages=sample_dict, # changed to lineages
                                                   double_step_validation=False)
    steps = steps * generation_time
    # some mig stuff
    mm = [x for x in ddb.demography.events if "MassMigration" in str(type(x))]
    census_size = ddb.population_size_trajectory(steps=steps)
    for m in reversed(mm):
        if m.proportion == 1.0:
            n = (steps > m.time)
            census_size[n, m.source] = census_size[n, m.dest]
        else:
            print(
                "Error: census size estimate requires that MassMigration proportion == 1.0")
            sys.exit(1)
    return coal_rate, P, steps, census_size


def gather_coal_rate(outfile, model, pops_dict, generation_time, steps):
    header = "method,population,DFE,annotations,year,Ne"
    with open(outfile[0], 'w') as coal:
        coal.write(f"{header}\n")
        for pop_id, pop in enumerate(pops_dict.keys()):
            coal_rate, P, step_rate, census_size = popn_coal_rate(model, pop, pops_dict, generation_time, steps)
            for i in range(len(step_rate)):
                coal.write(f"coal,{pop},none,none,{step_rate[i]},{1/(2*coal_rate[i])}\n")
                coal.write(f"census,{pop},none,none,{step_rate[i]},{census_size[i][pop_id]}\n")


def plot_compound_Ne_t(infile, outfile, ref_line="coal", colorby="annotations", log=True):
    """
    figure of N(t) for multiple methods
    set msmc w/ style="n_samp"
    """
    # load df
    df = pd.read_csv(infile, sep=",")
    df_ddb = pd.read_csv(Path(infile).parent.parent / "coal_estimated_Ne_t.csv", sep=",")
    df_ddb = df_ddb[df_ddb["method"] == ref_line]
    # plot params
    pop_order = np.sort(df['population'].unique())
    annot_order = np.sort(df[colorby].unique())[::-1]
    pal_dict = {pop:COLOURS[-(i+1)] for i, pop in enumerate(annot_order)}
    dfe_order = np.sort(df["DFE"].unique())[::-1]
    g = sns.relplot(data=df, x="year", y="Ne", col="population", row="DFE",
                    col_order=list(pop_order), row_order=list(dfe_order),
                    hue=colorby, hue_order=annot_order,
                    units="seed", kind="line", palette=pal_dict, 
                    height=3, aspect=2, estimator=None, errorbar="se", err_style="band",
                    facet_kws=dict(sharex=True, sharey=True), drawstyle='steps-pre')
    # add ref_line to all plots x pop name
    for ax in g.axes.flat:
        pop_ax = ax.title.get_text().split()[-1]
        df_ddb_pop = df_ddb.query(f"population == '{pop_ax}'")
        ax.plot(df_ddb_pop["year"], df_ddb_pop["Ne"], color="black")
    g.despine()
    #g.set(xlim=(0, max(df_ddb["year"])))
    if log:
        g.set(xscale="log", yscale="log")
    g.set_xlabels("time (years ago)")
    g.set_ylabels("population size")
    plt.savefig(f"{outfile}", bbox_inches='tight')


def plot_all_ne_estimates(infile, outfile, ref_line="coal", colorby="method", styleby="annotations", log=True):
    df = pd.read_csv(infile, sep=",")
    df_ddb = pd.read_csv(Path(infile).parent / "coal_estimated_Ne_t.csv", sep=",")
    df_ddb = df_ddb[df_ddb["method"] == ref_line]
    pop_order = np.sort(df["population"].unique())
    method_order = np.sort(df[colorby].unique())
    pal_dict = {pop:COLOURS[i] for i, pop in enumerate(method_order)}
    annot_order = np.sort(df[styleby].unique())[::-1]
    dfe_order = np.sort(df["DFE"].unique())[::-1]
    g = sns.relplot(data=df, x="year", y="Ne", row="DFE", col="population",
                    col_order=list(pop_order), row_order=list(dfe_order),
                    hue=colorby, hue_order=method_order, style=styleby, style_order=annot_order,
                    kind="line", drawstyle='steps-pre',
                    palette=pal_dict, errorbar="se", err_style="band", 
                    height=3, aspect=2, facet_kws=dict(sharex=True, sharey=True))
    # add ref_line to all plots x pop name
    for ax in g.axes.flat:
        pop_ax = ax.title.get_text().split()[-1]
        df_ddb_pop = df_ddb.query(f"population == '{pop_ax}'")
        ax.plot(df_ddb_pop["year"], df_ddb_pop["Ne"], color="black")
    g.despine()
    #g.set(xlim=(0, max(df_ddb["year"])))
    if log:
        g.set(xscale="log", yscale="log")
    g.set_xlabels("time (years ago)")
    g.set_ylabels("population size")
    plt.savefig(f"{outfile}", bbox_inches='tight')


# NOTE: is this all DFE below here?
def plot_all_dfe_results(input, output, mu, seq_len, nonneu_prop, pop_names=None):
    """
    Description:
        Plots the comparison of dadi, polyDFE, DFE-alpha and Grapes results for a single model.

    Arguments:
        input list: Names of input files.
        output list: Names of output files.
        mu float: Mutation rate per site per generation.
        seq_len int: Genomic sequence length.
        nonneu_prop float: Proportion of non-neutral mutations.
        pop_names list: Names of populations.
    """
    dadi_bestfits = _read_bestfits(input[0], len(pop_names))
    polydfe_bestfits = _read_bestfits(input[1], len(pop_names))
    dfe_alpha_bestfits = _read_bestfits(input[2], len(pop_names))
    grapes_bestfits = _read_bestfits(input[3], len(pop_names))

    # SLiM assumes genotype fitnesses: 1, 1+sh, 1+s
    # dadi and polyDFE assumes genotype fitnesses: 1, 1+2sh, 1+2s
    # s_SLiM = 2s_dadi = 2s_polyDFE

    # E(s_dadi) = shape*scale/(2*Na)
    ## Na = theta_ns/(4*mu*L_ns)
    # E(s_dadi) = shape*scale/(theta_ns/2*mu*L_ns)
    # See https://dadi.readthedocs.io/en/latest/user-guide/dfe-inference/
    dadi_shapes = [dadi_bestfits[i]['shape']
                   for i in range(len(dadi_bestfits))]
    dadi_Es = [dadi_bestfits[i]['shape']*dadi_bestfits[i]['scale'] /
               (dadi_bestfits[i]['theta']/(4*mu*seq_len*nonneu_prop)) for i in range(len(dadi_bestfits))]

    ## theta_bar = 4*Ne*mu
    ## S_d = 4*Ne*E(s)
    # E(s_polyDFE) = mu*S_d/theta_bar
    # See Fig. 1 in doi: 10.1534/genetics.117.300323
    polydfe_shapes = [polydfe_bestfits[i]['b']
                      for i in range(len(polydfe_bestfits))]
    polydfe_Es = [2*abs(mu*polydfe_bestfits[i]['S_d']/polydfe_bestfits[i]
                        ['theta_bar']) for i in range(len(polydfe_bestfits))]

    # Per https://doi.org/10.1534/genetics.107.080663 under Model
    # DFE alpha assumes genotype fitnesses are 1, 1-s/2, 1-s, so like SLiM
    dfe_alpha_shapes = [dfe_alpha_bestfits[i]['b']
                        for i in range(len(dfe_alpha_bestfits))]
    dfe_alpha_Es = [abs(dfe_alpha_bestfits[i]['Es'])
                    for i in range(len(dfe_alpha_bestfits))]
    
    grapes_shapes = [grapes_bestfits[i]['shape']
                     for i in range(len(grapes_bestfits))]
    grapes_Ne = [grapes_bestfits[i]['theta'] / (4*mu) for i in range(len(grapes_bestfits))]
    grapes_Es = [grapes_bestfits[i]['Es']/grapes_Ne[i] for i in range(len(grapes_bestfits))]

    fig = plt.figure(figsize=(6, 6), dpi=300)

    if pop_names is None:
        xtick_label = ['pop'+str(i) for i in np.arange(len(dadi_shapes))]
    else:
        xtick_label = pop_names
    plt.subplot(4, 2, 1)
    _plot_boxplots(plt, 'dadi', dadi_shapes, 0.19,
                   xtick_label, 'Shape', [0.1, 0.25])
    plt.subplot(4, 2, 2)
    _plot_boxplots(plt, 'dadi', dadi_Es, 0.014, xtick_label, '|E(s)|')

    plt.subplot(4, 2, 3)
    _plot_boxplots(plt, 'polyDFE', polydfe_shapes, 0.19,
                   xtick_label, 'Shape', [0.1, 0.25])
    plt.subplot(4, 2, 4)
    _plot_boxplots(plt, 'polyDFE', polydfe_Es, 0.014, xtick_label, '|E(s)|')

    plt.subplot(4, 2, 5)
    _plot_boxplots(plt, 'DFE-alpha', dfe_alpha_shapes,
                   0.19, xtick_label, 'Shape')
    plt.subplot(4, 2, 6)
    _plot_boxplots(plt, 'DFE-alpha', dfe_alpha_Es,
                   0.014, xtick_label, '|E(s)|')

    plt.subplot(4, 2, 7)
    _plot_boxplots(plt, 'grapes', grapes_shapes, 0.19,
                   xtick_label, 'Shape', print_xlabel=True)
    plt.subplot(4, 2, 8)
    _plot_boxplots(plt, 'grapes', grapes_Es, 0.014,
                   xtick_label, '|E(s)|', print_xlabel=True)

    fig.tight_layout()
    plt.savefig(output[0])


def _read_bestfits(input_files, pop_num):
    """
    Description:
        Combines all input files into a single data frame,
        then separates the single data frame into several data frames based on population ids.

    Arguments:
        input_files list: A list containing input files.
        pop_num int: Number of populations.

    Returns:
        bestfits list: A list containing data frames with bestfit parameters for different populations.
    """
    pop_ids = ['pop'+str(i) for i in np.arange(pop_num)]
    df_from_each_file = (pd.read_csv(f, sep="\t") for f in input_files)
    concatenated_df = pd.concat(df_from_each_file, ignore_index=True)
    bestfits = [concatenated_df[concatenated_df['pop_id'] == p]
                for p in pop_ids]

    return bestfits


def _plot_boxplots(plt, title, data, tval, xtick_label, ylabel, ylim=None, print_xlabel=False):
    """
    Description:
        Plots box plots for bestfit results from populations.

    Arguments:
        plt matplotlib.pyplot: Interface to matplotlib.
        title str: Title of the boxplot.
        data list: Bestfit results to be plotted.
        tval float: True value for the parameter.
        xtick_label list: X-axis tick labels.
        ylabel str: Y-axis label.
        ylim list: Y-axis limits.
        print_xlabel bool: True, print xlabel; False, not print xlabel.
    """
    plt.ylabel(ylabel)
    if print_xlabel:
        plt.xlabel('Population')
    plt.boxplot(data)
    plt.plot(plt.xlim(), [tval, tval], c='red')
    plt.xticks(np.arange(len(data))+1, xtick_label)
    if ylim is not None:
        plt.ylim(ylim)
    plt.title(title)


def _plot_binplots():
    """
    """
    pass


def _plot_sfsplots():
    """
    """
    pass
