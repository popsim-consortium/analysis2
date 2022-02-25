"""
Utilities for working with stairwayplot.
"""
import pathlib
import logging
import subprocess
import tempfile
import shutil
import concurrent.futures

import tskit
import allel
import plots
import numpy as np
import pandas as pd
import tqdm


def write_stairway_sfs(sequence_length, num_samples, sfs, path):
    """
    Writes the SFS to stairway plot format to the specified file.
    """
    # write stairwayplot input
    with open(path, "w") as out:
        # order is name,n_samples,sequence_length,lowest_sfs_bin,highest_sfs_bin
        print(
            "msp", num_samples, int(sequence_length), 1,
            num_samples - 1, sep="\t", file=out)
        # TODO can we use numpy to do this more efficiently?
        for x in sfs:
            print(int(x), end="\t", file=out)
        print(file=out)


class StairwayPlotRunner(object):
    """
    Run stairway plots.
    """
    def __init__(self, workdir, stairway_dir, java_exe="java"):
        self.workdir = pathlib.Path(workdir)
        shutil.rmtree(self.workdir, ignore_errors=True)
        self.workdir.mkdir(parents=True)
        stairway_path = pathlib.Path(stairway_dir)
        self.classpath = "{}:{}".format(stairway_path, stairway_path / "swarmops.jar")
        self.java_exe = java_exe

    def ts_to_stairway(self, ts_path, num_bootstraps=1, mask_file=None):
        """
        Converts the specified tskit tree sequence to text files used by
        stairway plot.
        """
        derived_counts_all = [[] for _ in range(num_bootstraps + 1)]
        total_length = 0
        num_samples = 0
        for i, ts_p in enumerate(ts_path):
            ts = tskit.load(ts_p)
            total_length += ts.sequence_length
            num_samples = ts.num_samples
            haps = ts.genotype_matrix()

            # Mapping mutation type IDs to class of mutation (e.g., neutral, non-neutral)
            mut_types = {}
            for dfe in ts.metadata["stdpopsim"]["DFEs"]:
                for mt in dfe["mutation_types"]:
                    mid = mt["slim_mutation_type_id"]
                    if not mid in mut_types:
                        mut_types[mid] = "neutral" if mt["is_neutral"] else "non_neutral"

            site_class = np.empty(ts.num_sites, dtype=object)

            # Finding neutral positions
            neu_positions = []
            non_neu_positions = []
            for j, s in enumerate(ts.sites()):
                mt = []
                for m in s.mutations:
                    for md in m.metadata["mutation_list"]:
                        mt.append(md["mutation_type"])
                        if set(mut_types[md["mutation_type"]]) == set("neutral"):
                            neu_positions.append(m.site)
                        if set(mut_types[md["mutation_type"]]) == set("non_neutral"):
                            non_neu_positions.append(m.site)
                site_class[j] = mut_types[mt[0]] if len(mt) == 1 else "double_hit" 
            assert sum(site_class == None) == 0

            # All SNP locs
            snp_locs =  [int(x.site.position) for x in ts.variants()]
            snp_locs = [snp_locs[index] for index in neu_positions]
            snp_locs_non_neutral = [snp_locs[index] for index in non_neu_positions]

            SFSs = []

            # Make it work, get ride off neutral positions that overlap the mask region
            retain = np.full(len(neu_positions), True)
            retain_non_neu = np.full(len(non_neu_positions), True)
            if mask_file:
                mask_table = pd.read_csv(mask_file, sep="\t", header=None)
                chrom = ts_p.split("/")[-1].split(".")[0]
                chrom = chrom.split("sim_")[1]
                sub = mask_table[mask_table[0] == chrom]
                mask_ints = pd.IntervalIndex.from_arrays(sub[1], sub[2])
                tmp_bool_neu = [mask_ints.contains(x) for x in snp_locs]
                tmp_bool_n_neu = [mask_ints.contains(x) for x in snp_locs_non_neutral]
                for i in range(len(tmp_bool_neu)):
                    if sum(tmp_bool_neu[i]) > 0:
                        retain[i] = False
                for i in range(len(tmp_bool_n_neu)):
                    if sum(tmp_bool_n_neu[i]) > 0:
                        retain_non_neu[i] = False
                total_length -= np.sum(mask_ints.length)

            neu_positions = list(neu_positions)

            # Extract neutral positions haplotypes
            haps_neu = haps[neu_positions,:]
            haps_non_neu = haps[retain_non_neu,:]            

            # append unmasked neutral SFS
            SFSs.append(allel.sfs(allel.HaplotypeArray(haps_neu).count_alleles()[:, 1])[1:])

            allele_counts = allel.HaplotypeArray(haps_neu).count_alleles()
            # get masked allele counts and append SFS
            allele_counts = allel.HaplotypeArray(haps_neu[retain,:]).count_alleles()

            # append masked neutral SFS
            SFSs.append(allel.sfs(allele_counts[:, 1])[1:])
            # append masked non neutral SFS
            SFSs.append(allel.sfs(allel.HaplotypeArray(haps_non_neu).count_alleles()[:, 1])[1:])

            sfs_path = ts_p+".sfs.pdf"
            # plotting masked and unmasked neutral SFSs and non neutral SFSs
            plots.plot_sfs(SFSs, sfs_path)


            # Bootstrap allele counts
            derived_counts_all[0].extend(allele_counts[:, 1])
            for j in range(1, num_bootstraps + 1):
                nsites = np.shape(allele_counts)[0]
                bootset = np.random.choice(np.arange(0, nsites, 1), nsites, replace=True)
                bootac = allele_counts[bootset, :]
                der_bootac = bootac[:, 1]
                derived_counts_all[j].extend(der_bootac)
        # Get the SFS minus the 0 bin and write output
        stairway_files = []
        for l in range(len(derived_counts_all)):
            sfs = allel.sfs(derived_counts_all[l])[1:]
            filename = self.workdir / f"sfs_{l}_.txt"
            write_stairway_sfs(total_length, num_samples, sfs, filename)
            stairway_files.append(filename)

        return stairway_files


    def _run_theta_estimation(self, input_file):
        """
        Runs stairway plot on the specified file, resulting in the output being written
        to the file input_file + ".addTheta".
        """
        num_runs = 1
        dim_factor = 5000
        # print("This is the input file")
        # print(input_file)
        cmd = (
            f"{self.java_exe} -cp {self.classpath} Stairway_plot_theta_estimation02 "
            f"{input_file} {num_runs} {dim_factor}")
        logging.info("Running:" + cmd)
        subprocess.run(cmd, shell=True, check=True)

    def run_theta_estimation(self, max_workers=None, show_progress=False):
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(self._run_theta_estimation, infile)
                for infile in self.workdir.glob(f"sfs_*.txt")]
            with tqdm.tqdm(total=len(futures), disable=not show_progress) as progress:
                for future in concurrent.futures.as_completed(futures):
                    progress.update()

    def run_summary(self, output_file, mutation_rate, generation_time):
        """
        Runs stairway plot summary files in the work dir, writing the
        output to output_file and with the given parameters.
        """
        # First we need to create a temporary directory for the files, stairway expects
        # only these files to be present in the directory.
        with tempfile.TemporaryDirectory() as tmpdir:
            for infile in self.workdir.glob(f"*_.txt.addTheta"):
                shutil.copy(infile, tmpdir)
            cmd = (
                f"{self.java_exe} -cp {self.classpath} "
                f"Stairway_plot_output_summary_commandline {tmpdir} "
                f"{mutation_rate} {generation_time} {output_file}")
            logging.info("Running:" + cmd)
            subprocess.run(cmd, shell=True, check=True)
