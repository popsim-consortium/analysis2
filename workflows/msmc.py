"""
Utilities for working with msmc
"""
import subprocess
import tskit
import numpy as np
import os
import pandas as pd
import glob
from pathlib import Path

def prune_tree_sequence(tree_sequence_path, pop_name):
    """
    take in a tree sequence, and a number of samples
    less than the number of samples in the tree,
    then simplify the tree sequence on that subset.
    """
    ts = tskit.load(tree_sequence_path)
    pop_id = [p.id for p in ts.populations() if p.metadata.get("name") == pop_name]
    pop_nodes = ts.samples(population=pop_id)
    ts = ts.simplify(samples=pop_nodes)
    return ts


def getIntervalOverlap(intA, intB):
    """
    returns the overlap between two pandas intervals
    """
    return max(0, min(intA.right, intB.right) - max(intA.left, intB.left))


def write_msmc_file(path, output, pop_name, mask_intervals=None):
    """
    take one .trees file and write out
    path.multihep.txt which acts as a single input to msmc

    This seems hacky atm, but let's getting working then hash out
    the details
    """
    ts = prune_tree_sequence(path, pop_name)
    dirr = os.path.dirname(path)
    filen = os.path.basename(path)
    sep = filen.split(".")
    chrom = sep[0].split("_")[1] # fixing sim_ name
    sep.insert(0, pop_name)
    fi = open(output, "w")
    prev = 0
    if mask_intervals is not None:
        mask_intervals = mask_intervals
        mask_intervals = pd.IntervalIndex.from_arrays(mask_intervals[:,0], mask_intervals[:,1])
        for var in ts.variants():
            cur = int(var.site.position)
            if cur > prev and sum(mask_intervals.contains(cur)) == 0:
                cur_int = pd.Interval(left=prev, right=cur)
                masked_bit = np.sum([getIntervalOverlap(cur_int, x) for x in
                                        mask_intervals])
                span = cur_int.length - masked_bit 
                geno = ''.join(map(str, var.genotypes))
                fi.write(f"{chrom}\t{cur}\t{span}\t{geno}\n")
                prev = cur
    else:
        for var in ts.variants():
            cur = int(var.site.position)
            #print("here")
            if cur > prev:
                geno = ''.join(map(str, var.genotypes))
                fi.write(f"{chrom}\t{cur}\t{cur-prev}\t{geno}\n")
                prev = cur

    fi.close()
    return None


def run_msmc_estimate(input_files, output_file, msmc_exec_loc, total_samples, num_genomes, iterations=1, ncores=1):
    """
    This is to run the msmc command and get estimates,
    It then will convert the scales times and pop sizes
    into generation years and population sizes.

    The final estimates will be written to
    input_file.final.txt
    """
    # TODO: change here, to num_samples and drop loop, if get wildcards.samps sorted in n_t.smk
    assert max(num_genomes) < total_samples * 2
    for nsamps in num_genomes:
        subset = np.random.choice(range(total_samples*2), nsamps, replace=False)
        haplotypes = ",".join(map(str, sorted(subset)))
        cmd = (f"{msmc_exec_loc} -r 0.25 -I {haplotypes} -i {iterations} -o {output_file}{nsamps}.trees.multihep.txt -t {ncores} {input_files}")
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(e.output)
            subprocess.run(cmd, shell=True, check=True)


def convert_msmc_output(results_file, outfile, mutation_rate, generation_time):
    """
    This function converts the output from msmc into a csv the will be read in
    for plotting comparison.

    MSMC outputs times and rates scaled by the mutation rate per basepair
    per generation. First, scaled times are given in units of the
    per-generation mutation rate. This means that in order to
    convert scaled times to generations, divide them by the mutation
    rate. In humans, we used mu=1e-8 per basepair per generation.
    To convert generations into years, multiply by the generation time,
    for which we used 10 years.

    To get population sizes out of coalescence rates, first take the
    inverse of the coalescence rate, scaledPopSize = 1 / lambda00. Then divide
    this scaled population size by 2*mu
    """
    infiles = glob.glob(str(Path(results_file).parent) + "/*.trees.multihep.txt.final.txt")
    with open(outfile, 'w') as out_fp:
        out_fp.write("pop\tyear\tNe\tn_samp\n")
        for infile in infiles:
            infile_path = Path(infile)
            n_samp = infile_path.name.split(".")[0]
            pop_name = infile_path.parent.parts[-1]
            in_fp = open(infile, "r")
            in_fp.readline()
            for line in in_fp:
                result = line.split()
                time = float(result[1])
                time_generation = time / mutation_rate
                time_years = time_generation * generation_time
                lambda00 = float(result[3])
                scaled_pop_size = 1 / lambda00
                size = scaled_pop_size / (2*mutation_rate)
                out_fp.write(f"{pop_name}\t{time_years}\t{size}\t{n_samp}\n")
