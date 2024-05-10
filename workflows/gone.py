"""
Utilities for working with scm++
"""
import subprocess
import tskit
import os
import numpy as np
import msprime
import pandas as pd

def params(gone_code, params):
    """
    modifies the GONE params file
    """
    with open(gone_code+"/INPUT_PARAMETERS_FILE") as params_file: text=params_file.readlines()
    text[4]=text[4].replace("PHASE=2", f"PHASE={str(params['gone_phase'])}")
    text[7]=text[7].replace("NGEN=2000", f"NGEN={str(params['gone_num_gens'])}")
    text[8]=text[8].replace("NBIN=400", f"NBIN={str(params['gone_num_bins'])}")
    text[12]=text[12].replace("maxNSNP=50000", f"maxNSNP={str(params['gone_max_snps'])}")
    #text[15]=text[15].replace("threads=-99", f"threads={str(params['gone_threads'])}")
    with open(f"{gone_code}/INPUT_PARAMETERS_FILE_TMP", "w") as params_file:
        for line in text:
            params_file.write(line)

    cmd = f"chmod u+x {gone_code}/PROGRAMMES/*"
    subprocess.run(cmd, shell=True, check=True)
    cmd = f"touch .params_edited"
    subprocess.run(cmd, shell=True, check=True)


def copy(gone_code, outpath, seed, threads):
    cwd = os.getcwd()
    g_list = [cwd, gone_code, outpath]
    cmd =  "".join(["ln -sf {}/{}/INPUT_PARAMETERS_FILE_TMP {}/INPUT_PARAMETERS_FILE;\n\n".format(*g_list)])
    cmd += "".join(["cat {}/script_GONE.sh | sed 's/num={}/num={}/g' | sed 's/Output_Ne_{}/gone_estimated_Ne.txt/g' > {}/script_GONE.sh;\n\n"
                    .format(gone_code, r'\$RANDOM', seed, r'\$FILE', outpath)])
    cmd += "".join([f"mkdir {outpath}/PROGRAMMES;\n\n"])
    cmd += "".join(["ln -sf {}/{}/PROGRAMMES/MANAGE_CHROMOSOMES2 {}/PROGRAMMES/MANAGE_CHROMOSOMES2;\n\n".format(*g_list)])
    cmd += "".join(["ln -sf {}/{}/PROGRAMMES/LD_SNP_REAL3 {}/PROGRAMMES/LD_SNP_REAL3;\n\n".format(*g_list)])
    cmd += "".join(["ln -sf {}/{}/PROGRAMMES/SUMM_REP_CHROM3 {}/PROGRAMMES/SUMM_REP_CHROM3;\n\n".format(*g_list)])
    cmd += "".join(["ln -sf {}/{}/PROGRAMMES/GONE {}/PROGRAMMES/GONE;\n\n".format(*g_list)])
    cmd += "".join(["ln -sf {}/{}/PROGRAMMES/GONEaverage {}/PROGRAMMES/GONEaverage;\n\n".format(*g_list)])
    cmd += "".join(["cat {}/PROGRAMMES/GONEparallel.sh | sed 's/threads={}/threads={}/g' | sed 's/{} {}\"/g' > {}/PROGRAMMES/GONEparallel.sh;\n\n"
                    .format(gone_code, r'\$(getconf _NPROCESSORS_ONLN)', str(threads), r'options_for_GONE=""/options_for_GONE="-sd', seed, outpath)])
    cmd += "".join([f"chmod u+x {outpath}/PROGRAMMES/GONEparallel.sh;\n\n"])
    cmd += "".join([f"touch {outpath}/.scripts_copied"])
    # print(cmd)
    subprocess.run(cmd, shell=True, check=True)


def ts2plink(ts_path, ped_file, map_file, species, pop_name, genetic_map, chromID, mask_intervals):
    """
    converts ts to plink format
    masks are the intervals to exclude
    """
    if type(mask_intervals) is not list:
        mask_intervals = [mask_intervals]
    if genetic_map is not None:
        gm_chr = [genetic_map.get_chromosome_map(chrms) for chrms in chromID]
    else:
        gm_chr = [species.get_contig(chrms).recombination_map for chrms in chromID]
    snp_counter = 1
    genomat_list = []
    # add to map file for chroms
    with open(map_file, "w") as mapfile:
        for i, ts_p in enumerate(ts_path):
            ts = tskit.load(ts_p)
            pop_id = [p.id for p in ts.populations() if p.metadata.get("name") == pop_name]
            pop_nodes = ts.samples(population=pop_id)
            ts = ts.simplify(samples=pop_nodes)
            genomat = ts.genotype_matrix() + 1 # it takes alleles (1,2)
            if mask_intervals[i] is not None:
                positions = ts.sites_position
                for interval in mask_intervals[i]:
                    masked_idx = np.where(np.logical_and(positions>=interval[0], positions<interval[1]))[0]
                    # missing_data state for plink is 0
                    genomat[masked_idx, :] = 0
            genomat_list.append(genomat)

            # make .map
            positions = ts.tables.sites.position
            # NOTE: this is very slow
            for snp in range(ts.genotype_matrix().shape[0]):
                bp = int(positions[snp])
                cM = gm_chr[i].get_cumulative_mass(bp) * 100
                #outline = [str(i+1), "snp"+str(snp_counter), "0", str(bp)] # test with uniform rec.
                outline = [str(i+1), "snp"+str(snp_counter), str(cM), str(bp)]
                mapfile.write(" ".join(outline) + "\n")
                snp_counter += 1

    # build inds
    inds = {}
    for sample in ts.samples():
        node = ts.node(sample)
        indID = node.individual
        if indID not in inds:
            inds[indID] = [sample]
        else:
            inds[indID].append(sample)

    # write after all chrom are added to genomat_list
    genomat = np.concatenate(genomat_list)
    with open(ped_file, "w") as pedfile:
        for ind in inds:
            outline=["1", "IND"+str(ind), "0", "0", "1", "-9"]
            nodes=inds[ind]
            for snp in range(genomat.shape[0]):
                outline.extend([str(genomat[snp][nodes[0]]), str(genomat[snp][nodes[1]])])
            pedfile.write(" ".join(outline) + "\n")
