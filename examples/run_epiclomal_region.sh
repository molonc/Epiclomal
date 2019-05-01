#! /bin/bash

# REPLACE THESE PATHS WITH YOUR PATHS
export PATH=/home/mandronescu/.local/centos6/miniconda3/bin:/home/mandronescu/.local/centos6/anaconda3/bin:$PATH
export LD_LIBRARY_PATH=/gsc/software/linux-x86_64-centos6/gcc-5.2.0/lib64/:$LD_LIBRARY_PATH
export R_LIBS=/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.5:/home/mandronescu/.local/centos6/R:$R_LIBS
export PYTHON=/home/mandronescu/.local/centos6/miniconda3/bin/python

# run EpiclomalRegion sequentially a specified number of times
# Should be at least 100 for real results, 300 or 1000 would be better.

for RUN in {0..50}
do
    echo "==================================="
    echo "Running EpiclomalRegion RUN $RUN"
    echo "==================================="    
    $PYTHON ../epiclomal/epiclomal.py Region-GeMM \
        --out_dir epiclomal_region/$RUN/ \
        --repeat_id $RUN \
        --true_prevalences 0.33_0.33_0.34 \
        --config_file inputs/config1.yaml \
        --true_clusters_file synthetic/data/true_clone_membership.tsv.gz \
        --methylation_file synthetic/data/data_incomplete.tsv.gz \
        --initial_clusters_file simple_hclust/initial_inputs.tsv.gz

done


# Possible input parameters for EpiclomalRegion

# $ python epiclomal/epiclomal.py Region-GeMM -h
# usage: Epiclomal Region-GeMM [-h] [--K K] --config_file CONFIG_FILE
#                              --methylation_file METHYLATION_FILE
#                              [--copynumber_file COPYNUMBER_FILE]
#                              [--regions_file REGIONS_FILE]
#                              [--initial_clusters_file INITIAL_CLUSTERS_FILE]
#                              [--true_clusters_file TRUE_CLUSTERS_FILE]
#                              [--true_prevalences TRUE_PREVALENCES]
#                              [--repeat_id REPEAT_ID] [--bulk_file BULK_FILE]
#                              [--slsbulk_file SLSBULK_FILE]
#                              [--slsbulk_iterations SLSBULK_ITERATIONS]
#                              [--out_dir OUT_DIR] [--mu_has_k MU_HAS_K]
#                              [--convergence_tolerance CONVERGENCE_TOLERANCE]
#                              [--max_num_iters MAX_NUM_ITERS] [--seed SEED]
#                              [--labels_file LABELS_FILE]
#                              [--Bishop_model_selection BISHOP_MODEL_SELECTION]
#                              [--check_uncertainty CHECK_UNCERTAINTY]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   --K K                 Max number of clusters. If given, use this one instead
#                         of the number from the config file.
#   --config_file CONFIG_FILE
#                         Path to YAML format configuration file.
#   --methylation_file METHYLATION_FILE
#                         Path to methylation data input file.
#   --copynumber_file COPYNUMBER_FILE
#                         Path to copy number input file.
#   --regions_file REGIONS_FILE
#                         Path to regions input file.
#   --initial_clusters_file INITIAL_CLUSTERS_FILE
#                         Start from these clusters instead of random clusters.
#   --true_clusters_file TRUE_CLUSTERS_FILE
#                         Path to the true_clusters_file, if known. If given,
#                         params.yaml will contain the V-measure for this
#                         prediction.
#   --true_prevalences TRUE_PREVALENCES
#                         A string with the true prevalences for all the
#                         clusters, e.g. 0.33_0.33_0.34
#   --repeat_id REPEAT_ID
#                         A number >= 0. If there is a column with this number
#                         (excluding the first column and starting from 0) in
#                         the initial_clusters_file, use that column as initial
#                         clusters, else use random initialization.
#   --bulk_file BULK_FILE
#                         A file with 3 columns: locus, #methylated reads,
#                         #unmethylated reads. The beta prior will be
#                         initialized with these values
#   --slsbulk_file SLSBULK_FILE
#                         A file with 3 columns: locus, #methylated reads,
#                         #unmethylated reads. This will be used to perform an
#                         SLS search that optimizes a bulk satisfaction score
#                         based on perturbation of uncertain cells.
#   --slsbulk_iterations SLSBULK_ITERATIONS
#                         The number of iterations for the SLSbulk procedure.
#   --out_dir OUT_DIR     Path where output files will be written.
#   --mu_has_k MU_HAS_K   True or False depending on whether we want mu to
#                         depend on k or not
#   --convergence_tolerance CONVERGENCE_TOLERANCE
#   --max_num_iters MAX_NUM_ITERS
#   --seed SEED           Set random seed so results can be reproduced. By
#                         default a random seed is chosen.
#   --labels_file LABELS_FILE
#                         Path of file with initial labels to use.
#   --Bishop_model_selection BISHOP_MODEL_SELECTION
#                         True or False depending on whether we want to apply
#                         Corduneanu_Bishop model selection
#   --check_uncertainty CHECK_UNCERTAINTY
#                         True or False depending on whether we want to check
#                         whether the uncertainty is estimated correctly