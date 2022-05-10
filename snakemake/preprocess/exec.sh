#!/bin/bash
#SBATCH --account=def-cdesouza
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH -o ./test.out # STDOUT

source /home/emiliano/projects/def-cdesouza/Lab/Epiclomal/venv/ec-env/bin/activate
snakemake --configfile config_smallwood.yaml


