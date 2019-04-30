#! /bin/bash

# REPLACE THESE PATHS WITH YOUR PATHS
export PATH=/home/mandronescu/.local/centos6/miniconda3/bin:/home/mandronescu/.local/centos6/anaconda3/bin:$PATH
export LD_LIBRARY_PATH=/gsc/software/linux-x86_64-centos6/gcc-5.2.0/lib64/:$LD_LIBRARY_PATH
export R_LIBS=/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.5:/home/mandronescu/.local/centos6/R:$R_LIBS
export RSCRIPT=/gsc/software/linux-x86_64-centos6/R-3.5.0/bin/Rscript

$RSCRIPT ../scripts/hclust.R \
    --evaluate_clustering_software ../epiclomal/evaluate_clustering.py \
    --true_prevalences 0.33_0.33_0.34 \
    --index ch \
    --max_k 10 \
    --output_directory simple_hclust \
    --methylation_file synthetic/data/data_incomplete.tsv.gz \
    --regions_file synthetic/data/regions_file.tsv.gz \
    --true_clone_membership_file synthetic/data/true_clone_membership.tsv.gz



# usage: hclust.R [-h] [--output_directory OUTPUT_DIRECTORY]
#                 [--methylation_file METHYLATION_FILE]
#                 [--regions_file REGIONS_FILE] [--max_k MAX_K]
#                 [--true_clone_membership_file TRUE_CLONE_MEMBERSHIP_FILE]
#                 [--true_prevalences TRUE_PREVALENCES]
#                 [--evaluate_clustering_software EVALUATE_CLUSTERING_SOFTWARE]
#                 [--index INDEX] [--impute IMPUTE]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   --output_directory OUTPUT_DIRECTORY
#                         Path to the output directory
#   --methylation_file METHYLATION_FILE
#                         Path to methylation data
#   --regions_file REGIONS_FILE
#                         Path to region coordinates
#   --max_k MAX_K         maximum number of clusters to be considered when
#                         cutting the tree
#   --true_clone_membership_file TRUE_CLONE_MEMBERSHIP_FILE
#                         Path to true clone membership file
#   --true_prevalences TRUE_PREVALENCES
#                         The true prevalence, for example 0.33_0.33_0.34
#   --evaluate_clustering_software EVALUATE_CLUSTERING_SOFTWARE
#                         Path to the file evaluate_clustering.py
#   --index INDEX         Index to be used to choose number of clusters, default
#                         = ch
#   --impute IMPUTE       If it is 1, impute with the average per region/locus,
#                         if it is 0 do nothing.