[![build Status](https://github.com/popsim-consortium/analysis2/actions/workflows/dry-run.yml/badge.svg?branch=main)](https://github.com/popsim-consortium/analysis2/actions)

# Analysis 2
Analysis of inference methods on standard population models including selection.
Here's how to get going, if you'd like to run this analysis

# Get the Analysis repo
Now clone the analysis2 repo, and install its dependencies
```console
$ git clone git@github.com:popsim-consortium/analysis2.git
$ cd analysis2/
```

# Set up your python environment to run the analysis
We recommend you start by creating a new `conda` environment for the analysis. This can be done using the command below, which will
create a new `conda` env called `analysis2`. Currently the workflow is targeted to run on python 3.9

```console
$ conda env create -f environment.yml
$ conda activate analysis2
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

```yaml
configfile: "workflows/config/snakemake/tiny_config.yaml"
```

Having set the config file, and perhaps edited to your wishes,
you are now ready to try a dry run.

# Perform a dry run of the workflow
To make sure that things are working, next run a dry run of the complete
workflow

```console
$ snakemake -c1 all -n
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

```console
$ snakemake -c M all
```

One can run the `clean_` targets of the workflow similarly.

## Running only a portion of the workflow
Sometimes the user only wants to run a subsection of the workflow. 
This is possible using `Snakemake` with the `--snakefile` option
along with the component workflows we have included. For instance,
to just perform the simulation steps of the workflow using 10 CPUs 
the user can say

```console
$ snakemake -c 10 --snakefile workflows/simulation.snake
```

and only that part of the analysis pipeline will run. 
We currently have 3 sub-workflows: `simulations.snake`
which does the simulations, `n_t.snake` which performs 
N(t) type analyses (e.g. `msmc`), and `dfe.snake` which
houses the portion of the workflow that does estimation
of the DFE.


## Running the workflow on a cluster
Currently we have provided two example `Snakemake` profiles
that allow a user to run the `analysis2` workflow on quite
easily. These can be found in `workflows/config/snakemake/oregon_profile`
and `workflows/config/snakemake/arizona_profile`.

### Cluster profile files
Each of those directories contains a single file, `config.yaml`,
that lays out the cluster specific settings needed to launch jobs.

For instance the `oregon_profile/config.yaml`, which is meant
to run using a cluster with a `slurm` scheduler looks like this

```yaml
cluster:
        mkdir -p logs/{rule} &&
        sbatch
                --partition=kern,kerngpu
                --account=kernlab
                --cpus-per-task={threads}
                --mem={resources.mem_mb}
                --time={resources.time}
                --job-name=smk-{rule}-{wildcards}
                --output=logs/{rule}/{rule}-{wildcards}-%j.out
default-resources:
        - time=60
        - mem_mb=5000
        - threads=1
restart-times: 3
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 500
keep-going: True
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
use-conda: True
```

to adopt this to a new `slurm` cluster a user would have to:
- change the `partition` value to appropriately named partitions
- change the `account` name

### Run command
At the command line this eases things tremendously. We can lauch the entire
workflow simply with

```console
$ snakemake --profile workflows/config/snakemake/oregon_profile/
```

`Snakemake` then will take care of all of the communication with the cluster,
launching jobs, and monitoring them for completeness. 

### Further reading on profiles
There are a lot of great examples on how to set up profiles for running
workflows on various cluster architectures. One excellent resource
is this repository of publicly available profiles https://github.com/snakemake-profiles/doc

