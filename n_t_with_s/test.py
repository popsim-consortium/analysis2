
"""
Snakefile for running N_t analysis on stdpopsim.

Simply running `snakemake` will run all analysis 
defined by the arguments above.

Currently, for each rep,
This will run stairway plot, smc++, and msmc
on the data resulting from simulations
on all chromosomes included in chrm_list
for estimates of N_t (Ne through time).
"""

import pathlib
import sys
import os
import numpy as np
import stdpopsim
import stairway
import smc
import msmc
import plots
import tskit
import dadi
from matplotlib import pyplot as plt

# ###############################################################################
# KNOBS - 
# ###############################################################################


# A seed to replicate results
# TODO mutation rates

config = {
    "seed" : 12345,
    "population_id" : 0,
    "num_samples_per_population" : [20, 0, 0],
    "num_sampled_genomes_msmc" : "2,8",
    "num_sampled_genomes_per_replicate": "2,8",
    "num_msmc_iterations" : 20,
    "replicates" : 1, 
    "species" : "HomSap",
    "model" : "OutOfAfricaArchaicAdmixture_5R19",
    "genetic_map" : "HapMapII_GRCh37",
    "chrm_list" : "chr22",
    "mask_file" : "masks/HapmapII_GRCh37.mask.bed",
    "dfe_id": "Gamma_K17",
    "selection":"all", 
    "annotation_id": "ensembl_havana_104_exons"
}


seeds = np.random.seed(config["seed"])

# This is the number of samples to simulate for within each population
# for each replicate
# TODO double check this is up to date with stdpopsim backend 
num_samples_per_population = config["num_samples_per_population"]
population_id = config["population_id"]
# Here is a list of sample sizes to run msmc on. 
# Each element counts as its own analysis
# so there will be "replicates" runs for each size
num_sampled_genomes_msmc = [int(x) for x in config["num_sampled_genomes_msmc"].split(",")]

# The number of msmc Baumwelch(?) iterations to run,
# typically 20
num_msmc_iterations = config["num_msmc_iterations"]

num_sampled_genomes_per_replicate = config["num_sampled_genomes_per_replicate"]

# The number of replicates of each analysis you would like to run
# For now leaving it a 1 just to get results quickly
replicates = config["replicates"]

# Define where selection will be imposed (whole genome or exons)
selection_list = [sel for sel in config["selection"].split(",")]
selection = selection_list[0]

# Where you would like all output files from analysis to live
output_dir = os.getcwd()

# The analysis species
species = stdpopsim.get_species(config["species"])

# The specific model you would like to run
if(config["model"] == "PiecewiseConstantSize"):
    model = stdpopsim.PiecewiseConstantSize(species.population_size)
    generation_time = species.generation_time
else:
    model = species.get_demographic_model(config["model"])
    # For plotting
    generation_time = model.generation_time

# The genetic map you would like to use.
# if value None is given default_recombination_rates are
# used with a flat map
genetic_map_id = config.get("genetic_map", None)

# The DFE id used for selection analyses
dfe_id = config.get("dfe_id", None)
# The GFF annotation id
annotation_id = config.get("annotation_id", None)

# This grabs the default mr from the first chromosome,
# Ultimitely This needs to be replaced with the weighted average
# of all chromosomes: This should be done in stdpopsim. 
mutation_rate = species.genome.mean_mutation_rate

chrms = config["chrm_list"]

# ###############################################################################
# GENERAL RULES & GLOBALS
# ###############################################################################


seed_array = np.random.random_integers(1,2**31,replicates)
genetic_map_downloaded_flag= ".genetic_map_downloaded"
msmc_exec = "../ext/msmc/build/msmc"
stairwayplot_code = "stairwayplot/swarmops.jar"

try:
    mask_file = config["mask_file"]
except KeyError:
    mask_file = None

contig = species.get_contig(chrms, genetic_map=genetic_map_id)
samples = model.get_samples(*num_samples_per_population)
engine = stdpopsim.get_engine("slim") ## here we use the slim engine
####### Here we are adding a given DFE
####### Adding DFE/selection to the simulation  
dfe = species.get_dfe(dfe_id)
selection = config["selection"]


if selection == "exons":
    ## Adding annotation only seletion on exon region
    exons = species.get_annotations(annotation_id)
    exon_intervals = exons.get_chromosome_annotations(chrms)
    contig.add_dfe(intervals=exon_intervals, DFE=dfe)

if selection == "all":
    # Adding selection to the whole contig
    #contig.add_dfe(intervals=np.array([[0, int(contig.length)]]), DFE=dfe)
    contig.add_dfe(intervals=np.array([[0, int(contig.length)]]), DFE=dfe)

if selection == "none":
    # if selection is not specified then a neutral model is expected
    contig = species.get_contig(chrms, genetic_map=genetic_map_id)

###### CHANGE SLIM SCALING FACTOR !!!!!!!!!!
# ts = engine.simulate(model, contig, samples, seed=seeds, slim_scaling_factor=10,
# slim_burn_in=10, slim_script=False)
# ts.dump(f"{chrms}_{selection}.trees")


def generate_fs_from_ts(ts, chrms, selection, sample_sets=None):

    def allele_counts(ts, sample_sets=None):
        if sample_sets is None:
            sample_sets = [ts.samples()]
        def f(x):
            return x
        return ts.sample_count_stat(sample_sets, f, len(sample_sets),
                                    span_normalise=False, windows='sites',
                                    polarised=True, mode='site', strict=False)

    # Mapping mutation type IDs to class of mutation (e.g., neutral, non-neutral)
    mut_types = {}
    for dfe in ts.metadata["stdpopsim"]["DFEs"]:
        for mt in dfe["mutation_types"]:
            mid = mt["slim_mutation_type_id"]
            if not mid in mut_types:
                mut_types[mid] = "neutral" if mt["is_neutral"] else "non_neutral"

    print("These are the mutation types", flush=True)
    print(mut_types, flush=True)

    site_class = np.empty(ts.num_sites, dtype=object)

    for j, s in enumerate(ts.sites()):
        mt = []
        for m in s.mutations:
            for md in m.metadata["mutation_list"]:
                mt.append(md["mutation_type"])
        site_class[j] = mut_types[mt[0]] if len(mt) == 1 else "double_hit" 
    assert sum(site_class == None) == 0
    
    print("Number of sites per class is:", flush=True)
    unique, counts = np.unique(site_class, return_counts=True)
    
    print(dict(zip(unique, counts)))
    freqs = allele_counts(ts, [sample_sets])
    freqs = freqs.flatten().astype(int)
    mut_classes = set(mut_types.values())
    mut_afs = {}
    # Feeding a dictionary with afs for each mutation type
    for mc in mut_classes:
        mut_afs[mc] = np.bincount(freqs[site_class == mc], minlength=len(sample_sets) + 1)

    print(mut_afs)
    # Extracting SFSs by summing each mutation type belonging two each category
    nonneu_fs = dadi.Spectrum(mut_afs["non_neutral"])
    neu_fs = dadi.Spectrum(mut_afs["neutral"])

    return [mut_afs["non_neutral"], mut_afs["neutral"]]

def plot_sfs(s, outfile):
    """
    Plot the SFS for this simulation
    """
    bins = [n + 1 for n in range(len(s[0]))]
    vals = []
    for i in range(len(s)):
        vals.append([int(x) for x in s[i]])
    if len(s) == 2:
        f, ax = plt.subplots(1, 2, sharey=True, tight_layout=True, figsize=(8, 3))
        ax[0].bar(bins, vals[0])
        ax[1].bar(bins, vals[1])
        ax[0].set_title("Non-neutral sites")
        ax[1].set_title("Neutral sites")
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


##### Cmputing the SFSs like Xin
ts = tskit.load(f"{chrms}_{selection}.trees")
sample_sets = ts.samples()

mut_afs = generate_fs_from_ts(ts, sample_sets=sample_sets, chrms=chrms, selection=selection)

print("These are the final results")
print(mut_afs)
plot_sfs(mut_afs, f"SFSs_{chrms}_{selection}.pdf")
