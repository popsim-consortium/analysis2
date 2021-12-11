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
            for md in m.metadata["mutation_list"]:
                mt.append(md["mutation_type"])
        if len(set(mt)) > 1:
            mut_type[j] = 3
        else:
            mut_type[j] = mt[0]

    freqs = allele_counts(ts, sample_sets)
    freqs = freqs.flatten().astype(int)
    mut_afs = np.zeros((len(sample_sets)+1, 3), dtype='int64')
    for k in range(3):
        mut_afs[:, k] = np.bincount(freqs[mut_type == k+1], minlength=len(sample_sets) + 1)

    return mut_afs

def generate_fs(ts, sample_sets, seq_len, neu_prop, nonneu_prop, output, format):
    """
    """
    mut_afs = _generate_fs_from_ts(ts, sample_sets)
    neu_fs = mut_afs[:,1]
    nonneu_fs = mut_afs[:,0]

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
    pass
def _generate_anavar_fs(ts, sample_sets, output):
    """
    """
    pass
