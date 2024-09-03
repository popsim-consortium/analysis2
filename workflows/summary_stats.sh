#config=workflows/config/snakemake/production_config_HomSap.yml
config=workflows/config/snakemake/production_config_PhoSin.yml

for i in {1..20};
do
    snakemake \
    --snakefile workflows/summary_stats.snake \
    --profile workflows/config/snakemake/oregon_profile/ \
    --configfile $config \
    --batch all=$i/20 \
    --groups run_diploshic_fvs=group0
done
