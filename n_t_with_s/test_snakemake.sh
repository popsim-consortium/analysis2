#!/bin/bash
#qsub -cwd -V -N simulation -l highp,time=100:00:00,h_data=16G -m bea test_snakemake.sh
/u/local/Modules/default/init/modules.sh
. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
conda activate popsim_env_test

snakemake -j 40 --config config="/u/home/m/mica20/project-kirk-bigdata/stdpopsim_analyses/demography_results/OutOfAfricaArchaicAdmixture_5R19"