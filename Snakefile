"""
Snakefile for running analysis 2 on stdpopsim.

Simply running snakemake -c all within this directory will run all analysis
defined by the modules below.

To run a single module (e.g. n_t workflow) one can do
snakemake -c 1 --snakefile workflows/n_t.snake --dry-run

Parameters are defined by the config.yaml file in 
workflows/config/snakemake/
"""

configfile: "workflows/config/snakemake/tiny_config.yaml"


module simulation_workflow:
    snakefile:
        "workflows/simulation.snake"
    config: config
use rule * from simulation_workflow as simulation_*

module dfe_workflow:
    snakefile:
        "workflows/dfe.snake"
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
rule clean_ext:
    shell:
        """
        rm -rf ext/GONE/ ext/grapes/ ext/stairwayplot/ ext/msmc/ 
        rm -rf ext/stairwayplot
        rm -rf ext/polyDFE
        """
rule clean_output:
    shell:
        """
        rm -rf {output_dir}
        """
rule clean_all:
    input:
        rules.clean_ext.input,
        rules.clean_output.input,
        #rules
    message: "Cleaning all"
    shell:
        """
        rm -rf ext/GONE/ ext/grapes/ ext/stairwayplot/ ext/msmc/ 
        rm -rf ext/stairwayplot
        rm -rf ext/polyDFE
        rm -rf {output_dir}
        """
    