# Analysis 2
Analysis of inference methods on standard population models including selection.
Here's how to get going, if you'd like to run this analysis

# Get the Analysis repo
Now clone the analysis2 repo, and install its dependencies
```
git clone https://github.com/popgensims/analysis2.git
cd analysis2/
```

# Set up your python environment to run the analysis
We recommend you start by creating a new `conda` environment for the analysis. 

```
conda create -n analysis2 python=3.9 \
    -c conda-forge --file conda-requirements.txt --yes
conda activate analysis2
pip install -r requirements.txt
```

For using `msmc` we need to download and compile it to play nice
with the conda environment that we have set up.
```
cd ext
git clone https://github.com/stschiff/msmc.git
cat msmc_makefile_stdpopsim_patch > msmc/Makefile && cd msmc && make
cd ../../
```


Further instructions can be currently found in each task directory.
A small example to simulate genomes with selection is described below.

------------
# Simulating selection with stdpopsim

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