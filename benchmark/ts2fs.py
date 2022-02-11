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

    mut_type = np.zeros(ts.num_sites)
    for j, s in enumerate(ts.sites()):
        mt = []
        for m in s.mutations:
            print(m)
            for md in m.metadata["mutation_list"]:
                mt.append(md["mutation_type"])
        if len(set(mt)) > 1:
            mut_type[j] = 3
        else:
            mut_type[j] = mt[0]

    freqs = allele_counts(ts, [sample_sets])
    freqs = freqs.flatten().astype(int)
    mut_afs = np.zeros((len(sample_sets)+1, 3), dtype='int64')
    for k in range(3):
        #mut_afs[:, k] = np.bincount(freqs[mut_type == k+1], minlength=len(sample_sets) + 1)
        mut_afs[:, k] = np.bincount(freqs[mut_type == k], minlength=len(sample_sets) + 1)

    return mut_afs

def generate_fs(ts, sample_sets, seq_len, neu_prop, nonneu_prop, output, format):
    """
    """
    mut_afs = _generate_fs_from_ts(ts, sample_sets)
    neu_fs = mut_afs[:,0]
    nonneu_fs = mut_afs[:,1]

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

if __name__ == '__main__':
    import stdpopsim
    import tskit
    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig('chr21')
    seq_len = contig.recombination_map.sequence_length
    neu_prop = 0.3
    nonneu_prop = 0.7
    seq_len = contig.recombination_map.sequence_length
    ts = tskit.load("/xdisk/rgutenk/xinhuang/projects/stdpopsim2/analysis2_dev/benchmark/results/simulated_data/878579711/sim_OutOfAfrica_3G09_HuberDFE_chr21.trees")
    samps = ts.samples()[:20]
    generate_fs(ts, samps, seq_len, neu_prop, nonneu_prop, 't', format='polydfe')

    #ts = tskit.load("results/simulated_data/OutOfAfrica_3G09/878579711/sim_OutOfAfrica_3G09_HuberDFE_chr21.trees")
    #samps = ts.samples()[:20]
    #generate_fs(ts, samps, seq_len, neu_prop, nonneu_prop, 's', format='polydfe')
