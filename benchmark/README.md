## Run specific tools

	snakemake --snakefile inference.snake -R --until combine_dadi_dfe_bestfits --profile ./config/snakemake/slurm/ --latency-wait 180
