"""
Utilities for working with scm++
"""
import subprocess
import tskit
import os
import numpy as np
import pandas as pd


def params(gone_code,params):
    """
    modifies the GONE params file
    """
    with open(gone_code+"/INPUT_PARAMETERS_FILE") as params_file: text=params_file.readlines()
    text[4]=text[4].replace("PHASE=2","PHASE="+str(params['gone_phase']))
    text[7]=text[7].replace("NGEN=2000","NGEN="+str(params['gone_num_gens']))
    text[8]=text[8].replace("NBIN=400","NBIN="+str(params['gone_num_bins']))
    text[12]=text[12].replace("maxNSNP=50000","maxNSNP="+str(params['gone_max_snps']))
    # text[15]=text[15].replace("threads=-99","threads="+str(params['gone_threads']))
    with open(gone_code+"/INPUT_PARAMETERS_FILE_TMP","w") as params_file: 
        for line in text:
            params_file.write(line)

    cmd = f"chmod u+x "+gone_code+"/PROGRAMMES/*"
    subprocess.run(cmd, shell=True, check=True)
    cmd = f"touch .params_edited"                  
    subprocess.run(cmd, shell=True, check=True)


def copy(gone_code,outpath,seed,threads):
    cwd = os.getcwd()
    cmd =  "".join(["ln -sf ",cwd,"/",gone_code,"/INPUT_PARAMETERS_FILE_TMP ",
                    outpath,"/INPUT_PARAMETERS_FILE;\n\n"])
    cmd += "".join(["cat " + gone_code + "/script_GONE.sh | sed s/\"num=\$RANDOM\"/\"RANDOM=", seed,"\\nnum=\$RANDOM\"/g > ", outpath, "/script_GONE.sh;\n\n"])
    cmd += "".join(["mkdir ", outpath, "/PROGRAMMES;\n\n"])
    cmd += "".join(["ln -sf ", cwd,"/",gone_code, "/PROGRAMMES/MANAGE_CHROMOSOMES2 ", outpath,"/PROGRAMMES/MANAGE_CHROMOSOMES2;\n\n"])
    cmd += "".join(["ln -sf ", cwd,"/",gone_code, "/PROGRAMMES/LD_SNP_REAL3 ", outpath,"/PROGRAMMES/LD_SNP_REAL3;\n\n"])
    cmd += "".join(["ln -sf ", cwd,"/",gone_code, "/PROGRAMMES/SUMM_REP_CHROM3 ", outpath,"/PROGRAMMES/SUMM_REP_CHROM3;\n\n"])
    cmd += "".join(["ln -sf ", cwd,"/",gone_code, "/PROGRAMMES/GONE ", outpath,"/PROGRAMMES/GONE;\n\n"])
    cmd += "".join(["ln -sf ", cwd,"/",gone_code, "/PROGRAMMES/GONEaverage ", outpath,"/PROGRAMMES/GONEaverage;\n\n"])
    cmd += "".join(["cat ", gone_code, "/PROGRAMMES/GONEparallel.sh | sed s/\"threads=\$(getconf _NPROCESSORS_ONLN)\"/\"threads=", str(threads) ,"\"/g | sed s/options_for_GONE=\\\"\\\"/options_for_GONE=\\\"-sd\ ", seed, "\\\"/g > ", outpath,"/PROGRAMMES/GONEparallel.sh;\n\n"])
    cmd += "".join(["chmod u+x ", outpath, "/PROGRAMMES/GONEparallel.sh;\n\n"])
    cmd += "".join(["touch ",outpath,"/.scripts_copied"])
    # print(cmd)
    subprocess.run(cmd, shell=True, check=True)

        
def ts2plink(ts, ped_file, map_file, gm_chr, chrID, mask_intervals=None):
    """
    converts ts to plink format
    masks are the intervals to exclude
    """
    ts = tskit.load(ts)
    inds = {}
    for sample in ts.samples():
        node = ts.node(sample)
        indID = node.individual
        if indID not in inds:
            inds[indID] = [sample]
        else:
            inds[indID].append(sample)

    # make .ped file
    genomat = ts.genotype_matrix() + 1 # it takes alleles (1,2)
    # mask if appropriate
    if mask_intervals is not None:
        positions = ts.sites_position
        for interval in mask_intervals:
            masked_idx = np.where(np.logical_and(positions>=interval[0], positions<interval[1]))[0]
            # missing_data state for plink is 0
            genomat[masked_idx, :] = 0

    with open(ped_file, "w") as pedfile:
        for ind in inds:
            outline=["1", "IND"+str(ind), "0", "0", "1", "-9"] 
            nodes=inds[ind]
            for snp in range(genomat.shape[0]):
                outline.extend([str(genomat[snp][nodes[0]]), str(genomat[snp][nodes[1]])])
            pedfile.write(" ".join(outline) + "\n")
    
    # make .map
    positions = ts.tables.sites.position
    with open(map_file, "w") as mapfile:
        for snp in range(genomat.shape[0]):
            bp = int(positions[snp])
            cM = gm_chr.get_cumulative_mass(bp)*100
            #outline = ["1", "snp"+str(snp), "0", str(bp)] # test with uniform rec.
            outline = ["1", "snp"+str(snp), str(cM), str(bp)]
            mapfile.write(" ".join(outline) + "\n")

