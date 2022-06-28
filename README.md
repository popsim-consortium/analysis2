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
create a new `conda` env called `analysis2`. Currently the workflow is targeted to run on python 3.8. We believe
there are issues when running python 3.9 or higher.

```
conda env create -f environment.yml
conda activate analysis2
```

For using `msmc` we need to download and compile it to play nice
with the conda environment that we have set up.
Also, we need DMD and GSL.
```
cd ext
    # GSL
wget https://ftp.gnu.org/gnu/gsl/gsl-2.7.tar.gz
tar -xzf gsl-2.7.tar.gz 
cd gsl-2.7
./configure --prefix=$PWD
make
make install
cd ../
    # DMD
wget https://s3.us-west-2.amazonaws.com/downloads.dlang.org/releases/2022/dmd.2.100.0.linux.tar.xz
tar -xf dmd.2.100.0.linux.tar.xz 
    # MSMC2
git clone https://github.com/stschiff/msmc2.git
gsl_path=$(echo $PWD/gsl-2.7/lib/ | sed s/"\/"/'\\\/'/g)
dmd_path=$(echo $PWD/dmd2/linux/bin64/dmd | sed s/"\/"/'\\\/'/g)
cd msmc2
mv Makefile Makefile_og
cat Makefile_og | sed s/dmd/$dmd_path/g | sed s/"\/usr\/local\/lib"/$gsl_path/g > Makefile
make
    #
cd ../../
```


For using `smc++` we need to download and build it.
```
cd ext
git clone https://github.com/popgenmethods/smcpp
cat smc_setup_stdpopsim_patch > smcpp/setup.py
cat smc_pyproject_stdpopsim_patch > smcpp/pyproject.toml
cd smcpp
pip install .
cd ../../
```

For using `DFE-alpha`, we need to download the program and extra data (5GB) from http://www.homepages.ed.ac.uk/pkeightl/dfe_alpha/download-dfe-alpha.html
```
cd ext
wget -c https://sourceforge.net/projects/dfe-alpha-k-e-w/files/dfe-alpha-release-2.16.tar.gz/download
mv download dfe-alpha-release-2.16.tar.gz
tar -xvf dfe-alpha-release-2.16.tar.gz
cat dfe_alpha_makefile_stdpopsim_patch > dfe-alpha-release-2.16/Makefile && cd dfe-alpha-release-2.16 && make
wget -c https://datashare.ed.ac.uk/bitstream/handle/10283/2730/data.tar.gz?sequence=1&isAllowed=y
tar -xvf data.tar.gz\?sequence\=1
cd ../../
```

For using `polyDFE` and `grapes`, we can download them from their GitHub repositories.
```
cd ext
git clone https://github.com/paula-tataru/polyDFE.git
cd polyDFE
chmod a+x polyDFE-2.0-linux-64-bit
wget -c https://github.com/BioPP/grapes/releases/download/v1.1.0/grapes-x86_64-bin-static-1.1.0-1.tar.gz
tar -xvf grapes-x86_64-bin-static-1.1.0-1.tar.gz
cd ..
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
