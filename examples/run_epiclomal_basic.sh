#! /bin/bash

# REPLACE THESE PATHS WITH YOUR PATHS
export PATH=/home/mandronescu/.local/centos6/miniconda3/bin:/home/mandronescu/.local/centos6/anaconda3/bin:$PATH
export LD_LIBRARY_PATH=/gsc/software/linux-x86_64-centos6/gcc-5.2.0/lib64/:$LD_LIBRARY_PATH
export R_LIBS=/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.5:/home/mandronescu/.local/centos6/R:$R_LIBS
export PYTHON=/home/mandronescu/.local/centos6/miniconda3/bin/python

# run EpiclomalBasic sequentially a specified number of times
# Should be at least 100 for real results, 300 or 1000 would be better.

for RUN in {0..50}
do
    echo "==================================="
    echo "Running EpiclomalBasic RUN $RUN"
    echo "==================================="    
    $PYTHON ../epiclomal/epiclomal.py Basic-GeMM \
        --out_dir epiclomal_basic/$RUN/ \
        --repeat_id $RUN \
        --true_prevalences 0.33_0.33_0.34 \
        --config_file inputs/config1.yaml \
        --true_clusters_file synthetic/data/true_clone_membership.tsv.gz \
        --methylation_file synthetic/data/data_incomplete.tsv.gz \
        --initial_clusters_file simple_hclust/initial_inputs.tsv.gz

done
