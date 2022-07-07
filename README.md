[![build Status](https://github.com/popsim-consortium/analysis2/actions/workflows/dry-run.yml/badge.svg?branch=main)](https://github.com/popsim-consortium/analysis2/actions)

# Analysis 2
Analysis of inference methods on standard population models including selection.
Here's how to get going, if you'd like to run this analysis

# Get the Analysis repo
Now clone the analysis2 repo, and install its dependencies
```
git clone git@github.com:popsim-consortium/analysis2.git
cd analysis2/
```

# Set up your python environment to run the analysis
We recommend you start by creating a new `conda` environment for the analysis. This can be done using the command below, which will
create a new `conda` env called `analysis2`. Currently the workflow is targeted to run on python 3.9

```
conda env create -f environment.yml
conda activate analysis2
```

# Set the config example
With the environment in place, the next step is to set the
workflow parameters using the a config file. 
`analysis2` currently ships with three example config
files, each found in `config/snakemake/`: `tiny_config.yaml`,
`config.yaml`, and `production-config.yaml`. Respectively
these represent a very small run, a small run, and the
final production settings used for the paper (TBD)

The workflow can be pointed at one of these config files
by editing the following line in the `Snakefile` file

```
configfile: "workflows/config/snakemake/tiny_config.yaml"
```

Having set the config file, and perhaps edited to your wishes,
you are now ready to try a dry run.

# Perform a dry run of the workflow
To make sure that things are working, next run a dry run of the complete
workflow

```
snakemake -c1 all -n
```

if the dry run checks out, you should be ready to run. 

# Run the workflow
You should now be set to run the complete workflow. This
will consist broadly of: 1) simulating the chromosomes
of interest, 2) downloading and installing the tools to be 
used in the anlaysis of the simulated data, 3) analysis of the
simulated tree sequences using the aforementioned tools, 
and 4) summarization of the analyses into figures. 

The Snakemake workflow has a number of targets:

- `all` -- run the complete analysis
- `clean_all` -- removes all simulations, downloads, and analysis
- `clean_ext` -- removes all downloaded external tools
- `clean_output` -- removes all simulation and analysis

To run the complete workflow on `M` cores use the following 

```
snakemake -c M all
```

One can run the `clean_` targets of the workflow similarly.

## Running only a portion of the workflow
Sometimes the user only wants to run a subsection of the workflow. 
This is possible using `Snakemake` with the `--snakefile` option
along with the component workflows we have included. For instance,
to just perform the simulation steps of the workflow using 10 CPUs 
the user can say

```
snakemake -c 30 --snakefile workflows/simulation.snake
```

and only that part of the analysis pipeline will run. 
We currently have 3 sub-workflows: `simulations.snake`
which does the simulations, `n_t.snake` which performs 
N(t) type analyses (e.g. `msmc`), and `dfe.snake` which
houses the portion of the workflow that does estimation
of the DFE.


------------------

A small example to simulate genomes with selection is described below.


# Simulating selection with stdpopsim

For each HomSap chromosome one can simulate selection in specific regions (e.g. coding region)
by using the function in stdpopsim and passed to
`stdpopsim.species.get_annotations` (see [the catalog](https://popsim-consortium.github.io/stdpopsim-docs/latest/tutorial.html#simulating-selection-on-exons) for details). Below is a complete example
that uses a given annotation to simulate selection on chr20.

```python
import stdpopsim

species = stdpopsim.get_species("HomSap")

# Adding a specific Distribution of fitness effects
dfe = species.get_dfe("Gamma_K17")

contig = species.get_contig("chr20")
model = species.get_demographic_model("OutOfAfrica_3G09")
samples = model.get_samples(100, 100, 100)  # YRI, CEU, CHB

# Extracting gene annottation
exons = species.get_annotations("ensembl_havana_104_exons")
exon_intervals = exons.get_chromosome_annotations("chr20")
contig.add_dfe(intervals=exon_intervals, DFE=dfe)

# Engine used for simulations
engine = stdpopsim.get_engine("slim")

# Simulate
ts = engine.simulate(
    model,
    contig,
    samples,
    seed=236,
    slim_scaling_factor=10,
    slim_burn_in=10,
)
```
