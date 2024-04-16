"""
Utilities for working with scm++
"""
import logging
import subprocess
import tskit
import numpy as np
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


def intervals2BedFile(mask_intervals, inter_mask_path, chr_name):
    with open(inter_mask_path, 'w') as bed:
        for int_s, int_e in mask_intervals:
            bed.write(f"{chr_name}\t{int_s - 1}\t{int_e}\n")


def write_smcpp_file(path, output, pop_name, num_sampled_genomes=2, mask_intervals=None):
    """
    Writes a smcpp input file given a treesequence
    """
    ts = prune_tree_sequence(path, pop_name)
    chr_name = Path(path).stem.split('_')[1]
    output_path = Path(output).parent
    outfile = f"{str(output_path)}/sim_{chr_name}.trees"
    vcf_file = f"{str(output_path)}/sim_{chr_name}.trees.vcf"
    mask_outfile = f"{str(output_path)}/sim_{chr_name}.trees.mask.bed"
    # write a vcf intermediate input
    with open(vcf_file, "w") as vcf:
        ts.write_vcf(vcf, contig_id=chr_name)  # site_mask=np.array(bool)
    # index/compress the vcf
    cmd = f"bgzip -f {vcf_file}"
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    vz_file = f"{vcf_file}.gz"
    cmd = f"tabix {vz_file}"
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    # write mask file
    if mask_intervals is not None:
        intervals2BedFile(mask_intervals, mask_outfile, chr_name)
        cmd = f"bgzip -f {mask_outfile}"
        logging.info("Running:" + cmd)
        subprocess.run(cmd, shell=True, check=True)
        cmd = f"tabix -p bed {mask_outfile}.gz"
        logging.info("Running:" + cmd)
        subprocess.run(cmd, shell=True, check=True)
        cmd_mask = f"smc++ vcf2smc --mask {mask_outfile}.gz -l {int(ts.sequence_length)} " 
    else:
        cmd_mask = f"smc++ vcf2smc -l {int(ts.sequence_length)} "
    # write command for smc++ vcf2smc
    dip_samps = ts.num_samples // 2
    inds = [f"tsk_{n}" for n in range(dip_samps)]
    inds = ','.join(inds)
    # NOTE: smc docs suggest pairs of haps for inference
    sampled_genomes = np.random.choice(range(dip_samps), (num_sampled_genomes-1, 2), replace=False)
    i=""  # dummy ext-name to fool snakemake
    for ind1, ind2 in sampled_genomes:
        cmd = cmd_mask + f"-d tsk_{ind1} tsk_{ind2} {vz_file} {outfile}{str(i).strip('0')}.smc.gz {chr_name} {pop_name}:{inds}"
        logging.info("Running:" + cmd)
        subprocess.run(cmd, shell=True, check=True)
        if i == "":
            i = 0.1
        else:
            i += 0.1


def run_smcpp_estimate(base, mutation_rate, ncores):
    """
    Runs smc++ estimate on the specified file, resulting in the output being written
    to the file input_file.final.jason".
    """
    cmd = (f"smc++ estimate --unfold --cores {ncores} {mutation_rate} *{base}")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    

def run_smcpp_plot(input_file, output_file, generation_time):
    """
    Runs smc++ plot on the specified file, resulting in the output being written
    to the file input_file.png".
    """
    cmd = (
        f"smc++ plot {input_file}.png {input_file} -g {generation_time} -c")
    logging.info("Running:" + cmd)
    subprocess.run(cmd, shell=True, check=True)
    cmd = (f"mv {input_file}.csv  {output_file}")
    subprocess.run(cmd, shell=True, check=True)
