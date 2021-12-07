Run `./get_intervals.py` to create interval files. This script:

 * downloads the HomSap GFF3 from ensembl,
 * extracts exons from the
   [Havana/Ensembl annotation merge](https://m.ensembl.org/info/genome/genebuild/annotation_merge.html),
 * and merges overlapping intervals (e.g. annotations on differing strands).

For each HomSap chromosome, a file is placed under `intervals/HomSap/`,
which can be loaded with `numpy.loadtxt(..., dtype=int)` and passed to
`stdpopsim.Contig.add_genomic_element_type()`. Below is a complete example
that uses these intervals to simulate background selection on chr22.

```python
import numpy as np
import stdpopsim

def KimDFE():
    """
    Return neutral and negative MutationType()s representing a human DFE.
    Kim et al. (2018), p.23, http://doi.org/10.1371/journal.pgen.1007741
    """
    neutral = stdpopsim.ext.MutationType()
    gamma_shape = 0.186  # shape
    gamma_mean = -0.01314833  # expected value
    h = 0.5  # dominance coefficient
    negative = stdpopsim.ext.MutationType(
        dominance_coeff=h,
        distribution_type="g",  # gamma distribution
        distribution_args=[gamma_mean, gamma_shape],
    )
    # neutral mutations have 0 proportion because they are not simulated by SLiM
    return {"mutation_types": [neutral, negative], "proportions": [0.0, 0.7]}

species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model("OutOfAfrica_3G09")
contig = species.get_contig("22")
contig.clear_genomic_mutation_types()
samples = model.get_samples(100, 100, 100)  # YRI, CEU, CHB

intervals = np.loadtxt("intervals/HomSap/ensembl_havana_exons_22.txt", dtype=int)
contig.add_genomic_element_type(intervals=intervals, **KimDFE())

# Simulate.
engine = stdpopsim.get_engine("slim")
ts = engine.simulate(
    model,
    contig,
    samples,
    seed=1234,
    slim_scaling_factor=10,
    slim_burn_in=10,
    # Set slim_script=True to print the script instead of running it.
    # slim_script=True,
)
```
