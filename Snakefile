
configfile: "workflows/config/snakemake/config.yaml"

module simulation_workflow:
    snakefile:
        "workflows/simulation.snake"
    config: config
use rule * from simulation_workflow as simulation_*

module dfe_workflow:
    snakefile:
        "workflows/inference.snake"
    config: config
use rule * from dfe_workflow as dfe_*

module n_t_workflow:
    snakefile:
        "workflows/n_t.snake"
    config: config
use rule * from n_t_workflow as n_t_*

output_dir = os.path.abspath(config["output_dir"])

# Define a new default target that collects all the targets from the modules.
rule all:
    input:
        rules.simulation_all.input,
        rules.dfe_all.input,
        rules.n_t_all.input,
       # rules
    default_target: True