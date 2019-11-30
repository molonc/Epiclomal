#! /bin/bash

# REPLACE THESE PATHS WITH YOUR PATHS
export PATH=/home/mandronescu/.local/centos6/miniconda3/bin:/home/mandronescu/.local/centos6/anaconda3/bin:$PATH
export LD_LIBRARY_PATH=/gsc/software/linux-x86_64-centos6/gcc-5.2.0/lib64/:$LD_LIBRARY_PATH
export R_LIBS=/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.5:/home/mandronescu/.local/centos6/R:$R_LIBS
export RSCRIPT=/gsc/software/linux-x86_64-centos6/R-3.5.0/bin/Rscript

# Evaluate EpiclomalBasic runs
$RSCRIPT ../scripts/eval_epiclomal.R \
    --regions_file synthetic/data/regions_file.tsv.gz \
    --output_dir results_basic \
    --input_dir epiclomal_basic/ \
    --true_clusters_file synthetic/data/true_clone_membership.tsv.gz \
    --true_epigenotypes_file synthetic/data/true_clone_epigenotypes.tsv.gz \
    --methylation_file synthetic/data/data_incomplete.tsv.gz \
    --model_name basic \
    --hdist_software ../scripts/hamming_distance.R \
    --visualization_software ../scripts/visualization.R



# $ Rscript scripts/eval_epiclomal.R -h
# usage: scripts/eval_epiclomal.R [-h] [--input_dir INPUT_DIR]
#                                 [--output_dir OUTPUT_DIR]
#                                 [--model_name MODEL_NAME]
#                                 [--hdist_software HDIST_SOFTWARE]
#                                 [--visualization_software VISUALIZATION_SOFTWARE]
#                                 [--methylation_file METHYLATION_FILE]
#                                 [--regions_file REGIONS_FILE]
#                                 [--true_clusters_file TRUE_CLUSTERS_FILE]
#                                 [--true_epigenotypes_file TRUE_EPIGENOTYPES_FILE]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   --input_dir INPUT_DIR
#                         Input directory of epiclomal results
#   --output_dir OUTPUT_DIR
#                         Output directory of the evaluation results
#   --model_name MODEL_NAME
#                         A name for the model
#   --hdist_software HDIST_SOFTWARE
#                         Full path of software that calculates hamming distance
#   --visualization_software VISUALIZATION_SOFTWARE
#                         Full path of software for visualization
#   --methylation_file METHYLATION_FILE
#                         Input data methylation file
#   --regions_file REGIONS_FILE
#                         Input regions file, has to be given even for basic
#   --true_clusters_file TRUE_CLUSTERS_FILE
#                         File with the true clusters, if known
#   --true_epigenotypes_file TRUE_EPIGENOTYPES_FILE
#                         File with the true epigenotypes, if known