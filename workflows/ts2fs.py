import dadi
import numpy as np
from numpy.random import default_rng
import tskit
import masks
import msprime

def mask_to_ratemap(annot_intervals, sequence_length, mutation_rate, proportion):
    """
    Convert an interval mask [a, b) to a mutation rate map

    Args:
    annot_intervals: np.array of shape (n_intervals, 2)
    sequence_length: int
    mutation_rate: float
    proportion: float -- proportion of mutations in interval that contribute to count
    """
    assert annot_intervals.shape[1] == 2
    bitmask = np.full(sequence_length, False)
    for (a, b) in annot_intervals:
        assert a >= 0, "Mask has negative coordinates"
        assert b <= sequence_length, "Mask exceeds sequence boundary"
        bitmask[int(a):int(b)] = True
    switch = np.logical_xor(bitmask[:-1], bitmask[1:])
    coords = np.append(0, np.flatnonzero(switch) + 1)
    values = bitmask[coords].astype(np.float64)
    coords = np.append(coords, sequence_length).astype(np.float64)
    values[values > 0] = mutation_rate * proportion
    return msprime.RateMap(position = coords, rate = values)

def get_expected_netural_subs(ts, mutation_rate):
    """
    Extract the expected number of neutral substitutions from a tree sequence

    Args:
    ts: tskit.TreeSequence
    mutation_rate: float
    """
    sequence_length = int(ts.sequence_length)
    exp_neutral_subs = 0.0
    for dfe in ts.metadata['stdpopsim']['DFEs']:
        for mtypes in dfe['mutation_types']:
            if mtypes['is_neutral']:
                annot_intervals = np.array(dfe['intervals'])
                if annot_intervals.shape[0] == 0:
                    annot_intervals = np.array([[0, sequence_length]])
                #assumes netural is first proportion in list
                proportion = dfe['proportions'][0]
                ratemap = mask_to_ratemap(annot_intervals, sequence_length, mutation_rate, proportion)
                slim_simulation_time = ts.metadata['SLiM']['tick']
                for t in ts.trees():
                    if t.num_edges > 0:
                        exp_neutral_subs += max(0, (slim_simulation_time - t.root) *
                                                (ratemap.get_cumulative_mass(t.interval.right) -
                                                 ratemap.get_cumulative_mass(t.interval.left)))
    return int(exp_neutral_subs)

def _generate_fs_from_ts(ts_dict, pop_idx, mask_file, annot, species, max_haplos):
    """
    Description:
        Converts tree sequences into 1d site frequency spectra.

    Arguments:
        ts_list dict: Dictionary of filepaths to tree sequences (values) for each chrm (key).
        pop_idx int: Population ID within the tree sequences.

    Returns:
        mut_afs numpy.ndarray: Frequency spectra for different types of mutations, summed over chrms.
    """
    def allele_counts(ts, sample_sets):
        """
        """
        def f(x):
            return x

        return ts.sample_count_stat(sample_sets, f, len(sample_sets),
                                    span_normalise=False, windows='sites',
                                    polarised=True, mode='site', strict=False)
    
    mut_afs = {"neutral":[], "non_neutral":[]}
    seq_len = 0
    neutral_anc_count = 0
    non_neutral_anc_count = 0
    neutral_prop = 0.3
    nonneu_prop = 0.7
    total_neutral_subs = 0
    for chrm in ts_dict:
        ts = tskit.load(ts_dict[chrm])
        contig = species.get_contig(chrm)
        chrom_length =  contig.recombination_map.sequence_length
        sample_sets = ts.samples(population=pop_idx)
        if max_haplos:
            sample_sets = sample_sets[:max_haplos]

        # apply masks
        if mask_file is not None:
            mask_intervals = masks.get_mask_from_file_dfe(mask_file, chrm)
            assert mask_intervals.shape[1] == 2
            mask_length = sum([interval[1] - interval[0] for interval in mask_intervals])
            ts = ts.delete_intervals(mask_intervals)

        # grab coding regions
        if (annot != "all_sites") and (annot != "none"):
            annotations = species.get_annotations(annot)
            annot_intervals = annotations.get_chromosome_annotations(chrm)
            ts = ts.keep_intervals(annot_intervals)
            exon_len = np.sum(annot_intervals[:,1]-annot_intervals[:,0])
            seq_len += exon_len
            non_neutral_anc_count += exon_len * nonneu_prop
            neutral_anc_count += exon_len * neutral_prop
        else:
            seq_len += chrom_length - mask_length
            non_neutral_anc_count += chrom_length * nonneu_prop
            neutral_anc_count += chrom_length * neutral_prop

        # Get count of substitutions for neutral SFS 
        mutation_rate = species.genome.mean_mutation_rate
        total_neutral_subs += get_expected_netural_subs(ts, mutation_rate)

        # Mapping mutation type IDs to class of mutation (e.g., neutral, non-neutral)
        mut_types = {}
        for dfe in ts.metadata["stdpopsim"]["DFEs"]:
            for mt in dfe["mutation_types"]:
                mids = mt["slim_mutation_type_id"]
                for mid in mids:
                    if not mid in mut_types:
                        mut_types[mid] = "neutral" if mt["is_neutral"] else "non_neutral"

        site_class = np.empty(ts.num_sites, dtype=object)
        for j, s in enumerate(ts.sites()):
            mt = []
            for m in s.mutations:
                for md in m.metadata["mutation_list"]:
                    mt.append(md["mutation_type"])
            site_class[j] = mut_types[mt[0]] if len(mt) == 1 else "double_hit"
        assert sum(site_class == None) == 0
        # print("Number of sites per class is:", flush=True)
        unique, counts = np.unique(site_class, return_counts=True)
        # print(dict(zip(unique, counts)), flush=True)

        # currently this is only working on single popn SFS
        freqs = allele_counts(ts, [sample_sets])
        freqs = freqs.flatten().astype(int)
        mut_classes = set(mut_types.values())

        # Feeding a dictionary with afs for each mutation type
        for mc in mut_classes:
            mut_afs[mc].append(np.bincount(freqs[site_class == mc], minlength=len(sample_sets) + 1))

    # sum over chrms
    mut_afs["neutral"] = np.array(mut_afs["neutral"])
    mut_afs["non_neutral"] = np.array(mut_afs["non_neutral"])
    mut_afs["neutral"] = np.sum(mut_afs["neutral"], axis=0)
    mut_afs["non_neutral"] = np.sum(mut_afs["non_neutral"], axis=0)

    #random draw of neutral subs from expected
    rng = default_rng(42)
    neutral_subs_draw = rng.poisson(total_neutral_subs)

    #add neutral subs to last bin 
    mut_afs["neutral"][-1] = neutral_subs_draw
    
    #add in ancestral counts
    mut_afs["neutral"][0] = neutral_anc_count
    mut_afs["non_neutral"][0] = non_neutral_anc_count
    
    # Remove mutated sites from first bin
    mut_afs["neutral"][0] = mut_afs["neutral"][0] - sum(mut_afs["neutral"][1:])
    mut_afs["non_neutral"][0] = mut_afs["non_neutral"][0] - sum(mut_afs["non_neutral"][1:])

    print(f"AFSs\nNeutral: {mut_afs['neutral']}\nNon-neutral: {mut_afs['non_neutral']}")
    return mut_afs,seq_len,len(sample_sets)


def generate_fs(ts_dict, pop_idx, mask_file, annot, species, output, format, max_haplos=None, seq_len=None, is_folded=False, **kwargs):
    """
    Description:
        Generates 1d site frequency spectra from tree sequences.

    Arguments:
        ts tskit.TreeSequence: A tree sequence.
        intervals np.ndarray: Coding regions used for generating frequency spectra.
        mask_intervals np.ndarray: Intervals removed from tree-sequences.
        sample_sets numpy.ndarray: A sample list.
        output list: Names of output files.
        format str: Format of output files. 
        is_folded bool: True, generates an unfolded frequency spectrum; False, generates a folded frequency spectrum.
        kwargs dict: Parameters for different DFE inference tools.
    """
    mut_afs,seq_len,sample_size = _generate_fs_from_ts(ts_dict, pop_idx, mask_file, annot, species, max_haplos)
    neu_fs = mut_afs["neutral"]
    nonneu_fs = mut_afs["non_neutral"]

    ## polyDFE v2.0 manual Page 1
    if is_folded and (format == 'polyDFE') : raise Exception('polyDFE only uses unfolded frequency spectra!') 
    if is_folded:
        neu_fs = _fold_fs(neu_fs)
        nonneu_fs = _fold_fs(nonneu_fs)

    if format == 'dadi': _generate_dadi_fs(neu_fs, nonneu_fs, output)
    elif format == 'polyDFE': _generate_polydfe_fs(neu_fs, nonneu_fs, output, seq_len, sample_size, **kwargs)
    elif format == 'DFE-alpha': _generate_dfe_alpha_fs(neu_fs, nonneu_fs, output, is_folded, **kwargs)
    elif format == 'grapes': _generate_grapes_fs(neu_fs, nonneu_fs, output, is_folded, seq_len, sample_size, **kwargs)

    else: raise Exception(f'{format} is not supported!')


def _fold_fs(fs):
    """
    Description:
        Converts a 1d unfolded frequency spectrum into a 1d folded frequency spectrum.

    Arguments:
        fs numpy.ndarray: 1d unfolded frequency spectrum.

    Returns:
        folded_fs numpy.ndarray: 1d folded frequency spectrum.
    """
    size = len(fs)
    folded_fs = [ fs[i] + fs[size-i-1] for i in range(int(size/2))]
    if size % 2 == 1: folded_fs += [ fs[int(size/2)] ]
    folded_fs = np.array(folded_fs)
  
    return folded_fs


def _generate_dadi_fs(neu_fs, nonneu_fs, output):
    """
    Description:
        Outputs frequency spectra for dadi.

    Arguments:
        neu_fs numpy.ndarray: Frequency spectrum for neutral mutations.
        nonneu_fs numpy.ndarray: Frequency spectrum for non-neutral mutations.
        output list: Names of output files.
    """
    neu_fs = dadi.Spectrum(neu_fs)
    nonneu_fs = dadi.Spectrum(nonneu_fs)

    neu_fs.to_file(output[0])
    nonneu_fs.to_file(output[1])


def _generate_polydfe_fs(neu_fs, nonneu_fs, output, seq_len, sample_size, **kwargs):
    """
    Description:
        Outputs frequency spectra for polyDFE.

    Arguments:
        neu_fs numpy.ndarray: Frequency spectrum for neutral mutations.
        nonneu_fs numpy.ndarray: Frequency spectrum for non-neutral mutations.
        output list: Names of output files.
        seq_len int: Length of the genomic sequence generates mutations.
        sample_size int: Number of haplotypes.
    Keyword Arguments:
        neu_prop float: Proportion of neutral mutations.
        nonneu_prop float: Proportion of non-neutral mutations
    """
    neu_len = round(seq_len * kwargs['neu_prop'])
    nonneu_len = round(seq_len * kwargs['nonneu_prop'])

    with open(output[0], 'w') as o:
        o.write(f"1 1 {sample_size}\n")
        o.write(" ".join([str(round(f)) for f in neu_fs[1:-1]]) + " " + str(neu_len) + "\n")
        o.write(" ".join([str(round(f)) for f in nonneu_fs[1:-1]]) + " " + str(nonneu_len) + "\n")


def _generate_dfe_alpha_fs(neu_fs, nonneu_fs, output, is_folded, **kwargs):
    """
    Description:
        Outputs frequency spectra for DFE-alpha.

    Arguments:
        neu_fs numpy.ndarray: Frequency spectrum for neutral mutations.
        nonneu_fs numpy.ndarray: Frequency spectrum for non-neutral mutations.
        output list: Names of output files.
        is_folded bool: True, generates a folded frequency spectrum; False, generates an unfolded frequency spectrum.
    
    Keyword Arguments (More details in DFE-alpha's manual):
        data_path_1 str: Directory containing the data files for one and two epoch models.
        data_path_2 str: Directory containing the data files for three epoch models.
        sfs_input_file str: Input file containing site frequency spectra.
        est_dfe_results_dir str: Directory containing results from est_dfe.
        est_dfe_demography_results_file str: Directory containing demography results from est_dfe.
        epochs int: 1, one-epoch model; 2, two-epoch model; 3, three-epoch model.
        site_class int: 0, analyse neutral SFS only; 1, analyse non-neutral SFS only; 2, analyse neutral and non-neutral SFS simultaneously; 2 is not allowed, if fold = 0 or epochs = 3.
        search_n2 int: Search for the best-fitting population size n2 in two-epoch model: 0, no; 1, yes.
    """
    def generate_config(data_path_1, data_path_2, sfs_input_file, est_dfe_results_dir, is_neutral, is_folded, est_dfe_demography_results_file=None, epochs=2):
        """
        """
        if is_folded: fold = 1
        else: fold = 0

        #if (site_class == 2) and (fold == 0): raise Exception('Cannot analyse neutral and non-neutral SFS simultaneously with unfolded SFS!')
        #if (site_class == 2) and (epochs == 3): raise Exception('Cannot analyse neutral and non-neutral SFS simultaneously with three-epoch model!')

        config = ''
        config += f'data_path_1\t{data_path_1}\n'
        config += f'data_path_2\t{data_path_2}\n'
        config += f'sfs_input_file\t{sfs_input_file}\n'
        if is_neutral:
            config += f'est_dfe_results_dir\t{est_dfe_results_dir}/neu\n'
            config += f'fold\t{fold}\n'
            config += f'epochs\t{epochs}\n'
            config += 'site_class\t0\n'
            config += 'search_n2\t1\n'
            config += 't2_variable\t1\n'
            config += 't2\t50\n'
        else:
            config += f'est_dfe_results_dir\t{est_dfe_results_dir}/nonneu\n'
            config += f'est_dfe_demography_results_file\t{est_dfe_demography_results_file}\n'
            config += f'fold\t{fold}\n'
            config += f'epochs\t{epochs}\n'
            config += 'site_class\t1\n'
            config += 'mean_s_variable\t1\n'
            config += 'mean_s\t-0.5\n'
            config += 'beta_variable\t1\n'
            config += 'beta\t0.5\n'
            config += 'p_additional\t0\n'
            config += 's_additional\t0\n'
        
        return config
  
    neu_config_out = output[0]
    nonneu_config_out = output[1]
    fs_out = output[2]

    neu_config = generate_config(kwargs['data_path_1'], kwargs['data_path_2'], kwargs['sfs_input_file'], kwargs['est_dfe_results_dir'], is_neutral=True, is_folded=is_folded)
    nonneu_config = generate_config(kwargs['data_path_1'], kwargs['data_path_2'], kwargs['sfs_input_file'], kwargs['est_dfe_results_dir'], is_neutral=False, is_folded=is_folded, est_dfe_demography_results_file=kwargs['est_dfe_demography_results_file'])

    ## DFE-alpha 2.16 manual Page 3
    ## If folded SFSs are provided by the user, then their length must be the number of alleles sampled +1,
    ## and the upper half of the SFS vectors should be zeros

    if is_folded:
       if len(neu_fs) % 2 == 1: upper_half = np.zeros(len(neu_fs))
       else: upper_half = np.zeros(len(neu_fs)-1)
       neu_fs = np.append(neu_fs, upper_half)

       if len(nonneu_fs) % 2 == 1: upper_half = np.zeros(len(nonneu_fs))
       else: upper_half = np.zeros(len(nonneu_fs)-1)
       nonneu_fs = np.append(nonneu_fs, upper_half)

    with open(neu_config_out, 'w') as o:
        o.write(neu_config)
    with open(nonneu_config_out, 'w') as o:
        o.write(nonneu_config)
    with open(fs_out, 'w') as o:
        o.write("1\n")
        o.write(str(len(neu_fs)-1)+"\n")
        o.write(" ".join([str(round(f)) for f in nonneu_fs]) + "\n")
        o.write(" ".join([str(round(f)) for f in neu_fs]) + "\n")

def _generate_grapes_fs(neu_fs, nonneu_fs, output, is_folded, seq_len, sample_size, **kwargs):
    """
    Description:
        Outputs frequency spectra for grapes.

    Arguments:
        neu_fs numpy.ndarray: Frequency spectrum for neutral mutations.
        nonneu_fs numpy.ndarray: Frequency spectrum for non-neutral mutations.
        output list: Names of output files.
        is_folded bool: True, generates a folded frequency spectrum; False, generates an unfolded frequency spectrum.
        seq_len int: Length of the genomic sequence generates mutations.
        sample_size int: Number of haplotypes.

    Keyword Arguments (Details in https://github.com/BioPP/grapes/blob/master/README.md):
        header str: A header.
        data_description str: Description for the data.
        neu_prop float: Proportion of neutral mutations.
        nonneu_prop float: Proportion of non-neutral mutations
    """
    neu_len = round(seq_len * kwargs['neu_prop'])
    neu_fs[0] = neu_len - sum(neu_fs[1:])
    nonneu_len = round(seq_len * kwargs['nonneu_prop'])
    nonneu_fs[0] = nonneu_len - sum(nonneu_fs[1:])

    with open(output[0], 'w') as o:
        o.write(kwargs['header']+"\n")
        if is_folded is not True: o.write("#unfolded\n")
        o.write(kwargs['data_description']+"\t")
        o.write(str(sample_size)+"\t")
        o.write(str(nonneu_len)+"\t")
        o.write("\t".join([str(f) for f in nonneu_fs[1:-1]])+"\t")
        o.write(str(neu_len)+"\t")
        o.write("\t".join([str(f) for f in neu_fs[1:-1]])+"\n")
