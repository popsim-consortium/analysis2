# Installing the requirements

```
conda create --name stdpopsim-analysis2 python=3.8
conda activate stdpopsim-analysis2
conda install -c conda-forge --file conda-requirements.txt
pip install -r requirements.txt
```

#Simulating selelection with stdpopsim

For each HomSap chromosome one can simulate selection in specific regions (e.g. coding region)
by using the function in stdpopsim and passed to
`stdpopsim.species.get_annotations` (see [the catalog](https://popsim-consortium.github.io/stdpopsim-docs/latest/tutorial.html#simulating-selection-on-exons) for details). Below is a complete example
that uses a given annotation to simulate selection on chr20.

```python
import numpy as np
import stdpopsim

species = stdpopsim.get_species("HomSap")

# Adding a specific Distribution of fitness effects
dfe = species.get_dfe("Gamma_K17")

contig = species.get_contig("chr20")
model = species.get_demographic_model("OutOfAfrica_3G09")
samples = model.get_samples(100, 100, 100)  # YRI, CEU, CHB

# Extracting gene annottation
exons = species.get_annotations("ensembl_havana_104_exons")
exon_intervals = exons.get_chromosome_annotations("chr20").astype("int")
contig.add_dfe(intervals=exon_intervals, DFE=dfe)

# Engine used for simulations
engine = stdpopsim.get_engine("slim")

# Simulate
ts = engine.simulate(
    model,
    contig,
    samples,
    seed=236,
    slim_scaling_factor=100,
    slim_burn_in=10,
)

```
