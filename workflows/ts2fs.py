import dadi
import numpy as np
import tskit


def _generate_fs_from_ts(ts, sample_sets=None):
    """
    Description:
        Converts tree sequences into 1d site frequency spectra.

    Arguments:
        ts tskit.TreeSequence: A tree sequence.
        sample_sets numpy.ndarray: A sample list.

    Returns:
        mut_afs numpy.ndarray: Frequency spectra for different types of mutations from the tree sequence.
    """
    if sample_sets is None:
        sample_sets = [ts.samples()]

    def allele_counts(ts, sample_sets):
        """
        """
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

    # print("These are the mutation types", flush=True)
    # print(mut_types, flush=True)

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
    mut_afs = {}
    for mc in mut_classes:
        mut_afs[mc] = np.bincount(freqs[site_class == mc], minlength=len(sample_sets) + 1)

    return mut_afs


def generate_fs(ts, sample_sets, output, format, coding_intervals=None, mask_intervals=None, is_folded=False, **kwargs):
    """
    Description:
        Generates 1d site frequency spectra from tree sequences.

    Arguments:
        ts tskit.TreeSequence: A tree sequence.
        intervals np.ndarray: Coding regsions used for generating frequency spectra.
        mask_intervals np.ndarray: Intervals removed from tree-sequences.
        sample_sets numpy.ndarray: A sample list.
        output list: Names of output files.
        format str: Format of output files. 
        is_folded bool: True, generates an unfolded frequency spectrum; False, generates a folded frequency spectrum.
        kwargs dict: Parameters for different DFE inference tools.
    """
    # If just a single sample set is provided, we need to wrap it
    # in a list to make it a list of sample setS
    #if not isinstance(sample_sets, list):
    #    sample_sets = [sample_sets]
    if coding_intervals is not None:
        ts = ts.keep_intervals(coding_intervals)
    if mask_intervals is not None:
        ts = ts.remove_intervals(mask_intervals)

    mut_afs = _generate_fs_from_ts(ts, sample_sets)
    neu_fs = mut_afs["neutral"]
    nonneu_fs = mut_afs["non_neutral"]

    ## polyDFE v2.0 manual Page 1
    if is_folded and (format == 'polyDFE') : raise Exception('polyDFE only uses unfolded frequency spectra!') 
    if is_folded:
        neu_fs = _fold_fs(neu_fs)
        nonneu_fs = _fold_fs(nonneu_fs)

    if format == 'dadi': _generate_dadi_fs(neu_fs, nonneu_fs, output)
    elif format == 'polyDFE': _generate_polydfe_fs(neu_fs, nonneu_fs, output, **kwargs)
    elif format == 'DFE-alpha': _generate_dfe_alpha_fs(neu_fs, nonneu_fs, output, is_folded, **kwargs)
    elif format == 'grapes': _generate_grapes_fs(neu_fs, nonneu_fs, output, is_folded, **kwargs)
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


def _generate_polydfe_fs(neu_fs, nonneu_fs, output, **kwargs):
    """
    Description:
        Outputs frequency spectra for polyDFE.

    Arguments:
        neu_fs numpy.ndarray: Frequency spectrum for neutral mutations.
        nonneu_fs numpy.ndarray: Frequency spectrum for non-neutral mutations.
        output list: Names of output files.
    
    Keyword Arguments:
        sample_size int: Number of haplotypes.
        seq_len int: Length of the genomic sequence generates mutations.
        neu_prop float: Proportion of neutral mutations.
        nonneu_prop float: Proportion of non-neutral mutations
    """
    neu_len = round(kwargs['seq_len'] * kwargs['neu_prop'])
    nonneu_len = round(kwargs['seq_len'] * kwargs['nonneu_prop'])

    with open(output[0], 'w') as o:
        o.write(f"1 1 {kwargs['sample_size']}\n")
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


def _generate_grapes_fs(neu_fs, nonneu_fs, output, is_folded, **kwargs):
    """
    Description:
        Outputs frequency spectra for grapes.

    Arguments:
        neu_fs numpy.ndarray: Frequency spectrum for neutral mutations.
        nonneu_fs numpy.ndarray: Frequency spectrum for non-neutral mutations.
        output list: Names of output files.
        is_folded bool: True, generates a folded frequency spectrum; False, generates an unfolded frequency spectrum.

    Keyword Arguments (Details in https://github.com/BioPP/grapes/blob/master/README.md):
        header str: A header.
        data_description str: Description for the data.
        sample_size int: Number of haplotypes.
        seq_len int: Length of the genomic sequence generates mutations.
        neu_prop float: Proportion of neutral mutations.
        nonneu_prop float: Proportion of non-neutral mutations
    """
    neu_len = round(kwargs['seq_len'] * kwargs['neu_prop'])
    nonneu_len = round(kwargs['seq_len'] * kwargs['nonneu_prop'])

    with open(output[0], 'w') as o:
        o.write(kwargs['header']+"\n")
        if is_folded is not True: o.write("#unfolded\n")
        o.write(kwargs['data_description']+"\t")
        o.write(str(kwargs['sample_size'])+"\t")
        o.write(str(nonneu_len)+"\t")
        o.write("\t".join([str(f) for f in nonneu_fs[1:-1]])+"\t")
        o.write(str(neu_len)+"\t")
        o.write("\t".join([str(f) for f in neu_fs[1:-1]])+"\n")
