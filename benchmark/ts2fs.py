import dadi
import numpy as np
import tskit


def _generate_fs_from_ts(ts, sample_sets=None):
    """
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

    freqs = allele_counts(ts, sample_sets)
    freqs = freqs.flatten().astype(int)
    mut_classes = set(mut_types.values())
    # Feeding a dictionary with afs for each mutation type
    mut_afs = {}
    for mc in mut_classes:
        mut_afs[mc] = np.bincount(freqs[site_class == mc], minlength=len(sample_sets) + 1)
    return mut_afs

def generate_fs(ts, sample_sets, seq_len, neu_prop, nonneu_prop, output, format):
    """
    """
    # If just a single sample set is provided, we need to wrap it
    # in a list to make it a list of sample setS
    if not isinstance(sample_sets, list):
        sample_sets = [sample_sets]
    mut_afs = _generate_fs_from_ts(ts, sample_sets)
    neu_fs = mut_afs["neutral"]
    nonneu_fs = mut_afs["non_neutral"]

    if format == 'dadi': _generate_dadi_fs(neu_fs, nonneu_fs, output)
    elif format == 'polydfe': _generate_polydfe_fs(neu_fs, nonneu_fs, seq_len, neu_prop, nonneu_prop, output)
    elif format == 'dfe-alpha': _generate_dfe_alpha_fs(ts, sample_sets, output)
    elif format == 'anavar': _generate_anavar_fs(ts, sample_sets, output)

def _generate_dadi_fs(neu_fs, nonneu_fs, output):
    """
    """
    neu_fs = dadi.Spectrum(neu_fs)
    nonneu_fs = dadi.Spectrum(nonneu_fs)
    neu_fs.to_file(output[0])
    nonneu_fs.to_file(output[1])

def _generate_polydfe_fs(neu_fs, nonneu_fs, seq_len, neu_prop, nonneu_prop, output):
    """
    """
    neu_len = round(seq_len * neu_prop)
    nonneu_len = round(seq_len * nonneu_prop)
    with open(output[0], 'w') as o:
        o.write("1 1 20\n")
        o.write(" ".join([str(round(f)) for f in neu_fs[1:-1]]) + " " + str(neu_len) + "\n")
        o.write(" ".join([str(round(f)) for f in nonneu_fs[1:-1]]) + " " + str(nonneu_len) + "\n")

def _generate_dfe_alpha_fs(ts, sample_sets, output):
    """
    """
    ## For DFE-alpha
#        dfe_alpha_neu_config = [
#            'data_path_1\t/xdisk/rgutenk/xinhuang/software/dfe-alpha-release-2.16/data',
#            'data_path_2\tdfe-alpha-three-epoch',
#            'sfs_input_file\t'+output_dir+"/dfe-alpha/"+wildcards.seeds+"/"+wildcards.chrms+"/pop"+wildcards.ids+"/pop"+wildcards.ids+".dfe-alpha.fs",
#            'est_dfe_results_dir\t'+output_dir+"/dfe-alpha/"+wildcards.seeds+"/"+wildcards.chrms+"/pop"+wildcards.ids+"/neu",
#            'fold\t0',
#            'epochs\t2',
#            'site_class\t0',
#            'search_n2\t1',
#            't2_variable\t1',
#            't2\t50'
#        ]
#        dfe_alpha_nonneu_config = [
#            'data_path_1\t/xdisk/rgutenk/xinhuang/software/dfe-alpha-release-2.16/data',
#            'data_path_2\tdfe-alpha-three-epoch',
#            'sfs_input_file\t'+output_dir+"/dfe-alpha/"+wildcards.seeds+"/"+wildcards.chrms+"/pop"+wildcards.ids+"/pop"+wildcards.ids+".dfe-alpha.fs",
#            'est_dfe_results_dir\t'+output_dir+"/dfe-alpha/"+wildcards.seeds+"/"+wildcards.chrms+"/pop"+wildcards.ids+"/nonneu",
#            'est_dfe_demography_results_file\t'+output_dir+"/dfe-alpha/"+wildcards.seeds+"/"+wildcards.chrms+"/pop"+wildcards.ids+"/neu/est_dfe.out",
#            'fold\t0',
#            'epochs\t2',
#            'site_class\t1',
#            'mean_s_variable\t1',
#            'mean_s\t-0.5',
#            'beta_variable\t1',
#            'beta\t0.5',
#            'p_additional\t0',
#            's_additional\t0',
#        ]
#        with open(output[0], 'w') as o:
#            o.write("\n".join(dfe_alpha_neu_config)+"\n")
#        with open(output[1], 'w') as o:
#            o.write("\n".join(dfe_alpha_nonneu_config)+"\n")
#        with open(output[2], 'w') as o:
#            o.write("1\n")
#            o.write(str(len(neu_fs)-1)+"\n")
#            o.write(" ".join([str(round(f)) for f in nonneu_fs.data]) + "\n")
#            o.write(" ".join([str(round(f)) for f in neu_fs.data]) + "\n")
    pass

def _generate_anavar_fs(ts, sample_sets, output):
    """
    """
        ## For anavar
        #anavar_alg_cmd = [
        #    '[algorithm_commands]',
        #    'search_algorithm: NLOPT_LD_SLSQP',
        #    'maxeval: 100000',
        #    'maxtime: 600',
        #    'num_searches: 100',
        #    'nnoimp: 2',
        #    'maximp: 3',
        #    'optional: false'
        #]
        #anavar_mdl_cmd = [
        #    '[model_commands]',
        #    'model: SNP_1',
        #    'n: ' + str(nonneu_fs.sample_sizes[0]),
        #    'm: ' + str(seq_len),
        #    'folded: false',
        #    'sfs: ' + ", ".join([str(round(f)) for f in nonneu_fs.data[1:-1]]),
        #    'dfe: continuous',
        #    'distribution: reflected_gamma',
        #    'theta_range: 0.1, 1000',
        #    'shape_range: 0.1, 100',
        #    'scale_range: 0.1, 100',
        #    'e_range: 0, 1',
        #    'constraint: no_pol_error',
        #    'optional: false'
        #]

        #with open(output[4], 'w') as o:
        #    o.write("\n".join(anavar_alg_cmd)+"\n\n\n")
        #    o.write("\n".join(anavar_mdl_cmd)+"\n")

        ## For prfreq
        #with open(output[5], 'w') as o:
        #    for d in nonneu_fs.data[1:]:
        #        o.write(f"{d}\n")

    pass
