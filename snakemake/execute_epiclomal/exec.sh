#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH -o ./test.out # STDOUT

source /home/emiliano/projects/def-cdesouza/Lab/Epiclomal/venv/ec-env/bin/activate
snakemake --configfile /home/emiliano/projects/def-cdesouza/Lab/Epiclomal/snakemake/execute_epiclomal/config_Smallwood2014.yaml
        