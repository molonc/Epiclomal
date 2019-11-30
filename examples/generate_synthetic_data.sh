#! /bin/bash

# REPLACE THESE PATHS WITH YOUR PATHS
export PATH=/home/mandronescu/.local/centos6/miniconda3/bin:/home/mandronescu/.local/centos6/anaconda3/bin:/gsc/software/linux-x86_64-centos6/R-3.5.0/bin/:$PATH
export LD_LIBRARY_PATH=/gsc/software/linux-x86_64-centos6/gcc-5.2.0/lib64/:$LD_LIBRARY_PATH
export R_LIBS=/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.5:/home/mandronescu/.local/centos6/R:/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.5:/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.5:$R_LIBS
export RSCRIPT=/gsc/software/linux-x86_64-centos6/R-3.5.0/bin/Rscript

# see below all the options of the generate_synthetic.R script
$RSCRIPT ../scripts/generate_synthetic.R \
    --saveall=0 \
    --verbose=0 \
    --clone_prevalence=0.33_0.33_0.34 \
    --read_size=10_2 \
    --region_size_type=multinomial_equal \
    --genotype_prob=dirichlet \
    --num_regions=10 \
    --num_clones=3 \
    --num_cells=100 \
    --plot_data=1 \
    --seed=31 \
    --prop_add_var=0_0.5 \
    --prop_cpg_flip=1 \
    --num_samples=1 \
    --missing_probability=0.6 \
    --error_probability=0.001_0.001 \
    --num_loci=1000 \
    --given_dir_complete=1 \
    --phylogenetic_generation=1 \
    --percent_regions_dirichlet_param=0.8 \
    --dirichlet_param_genotype_prob=90_10 \
    --output_dir=synthetic \
    --visualization_software ../scripts/visualization.R
    

# $ Rscript ../scripts/generate_synthetic.R -h
# usage: generate_synthetic.R [-h] [--read_size READ_SIZE]
#                             [--num_samples NUM_SAMPLES] [--num_loci NUM_LOCI]
#                             [--num_clones NUM_CLONES] [--num_cells NUM_CELLS]
#                             [--clone_prevalence CLONE_PREVALENCE]
#                             [--error_probability ERROR_PROBABILITY]
#                             [--missing_probability MISSING_PROBABILITY]
#                             [--genotype_prob GENOTYPE_PROB]
#                             [--dirichlet_param_genotype_prob DIRICHLET_PARAM_GENOTYPE_PROB]
#                             [--percent_regions_dirichlet_param PERCENT_REGIONS_DIRICHLET_PARAM]
#                             [--num_regions NUM_REGIONS]
#                             [--region_size_type REGION_SIZE_TYPE]
#                             [--output_dir OUTPUT_DIR]
#                             [--given_dir_complete GIVEN_DIR_COMPLETE]
#                             [--plot_data PLOT_DATA]
#                             [--visualization_software VISUALIZATION_SOFTWARE]
#                             [--prop_add_var PROP_ADD_VAR]
#                             [--bulk_depth BULK_DEPTH]
#                             [--phylogenetic_generation PHYLOGENETIC_GENERATION]
#                             [--prop_cpg_flip PROP_CPG_FLIP] [--seed SEED]
#                             [--verbose VERBOSE] [--saveall SAVEALL]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   --read_size READ_SIZE
#                         mean and std dev for the number of CpGs per read when
#                         generating data considering a read-based approach,
#                         default is no read based approach, that is 1_0
#   --num_samples NUM_SAMPLES
#                         Number different samples
#   --num_loci NUM_LOCI   Number of loci
#   --num_clones NUM_CLONES
#                         Number of clones
#   --num_cells NUM_CELLS
#                         Number of cells per sample, starting with sample 1,
#                         then sample 2, etc
#   --clone_prevalence CLONE_PREVALENCE
#                         Probability of clone prevalences for all samples
#   --error_probability ERROR_PROBABILITY
#                         Error probability for each methylation state
#   --missing_probability MISSING_PROBABILITY
#                         Missing data probability, it could be one probability
#                         for all cells or different ones
#   --genotype_prob GENOTYPE_PROB
#                         dirichlet (genotype probabilities are draws from a
#                         dirichlet distribution) or 0.5_fixed (all genotype
#                         probabilities fixed to 0.5)
#   --dirichlet_param_genotype_prob DIRICHLET_PARAM_GENOTYPE_PROB
#                         Dirichlet parameters to draw genotype probabilities
#                         for each region r and clone k
#   --percent_regions_dirichlet_param PERCENT_REGIONS_DIRICHLET_PARAM
#                         The percentage of the regions that have methylation
#                         drawn from the above distribution. For example if I
#                         want 90 percent of the regions to be hypermethylated
#                         and 10 percent to be hypomethylated, I will set
#                         dirichlet_param_genotype_prob=99_1 and
#                         percent_regions_dirichlet_param=0.9
#   --num_regions NUM_REGIONS
#                         Number of regions
#   --region_size_type REGION_SIZE_TYPE
#                         uniform, multinomial_equal or multinomial_nonequal.
#                         Fixed generated from uniform (from 1 to nloci),
#                         multinomial_equal (with prob 1/nregions) or
#                         multinomial_nonequal (currently this has hard-coded
#                         probabilities)
#   --output_dir OUTPUT_DIR
#                         Entire or just the beginning of the output directory
#                         file
#   --given_dir_complete GIVEN_DIR_COMPLETE
#                         If this is 0, it creates a long output dir name with
#                         the input parameters, if it is 1, the output dir is
#                         output_dir
#   --plot_data PLOT_DATA
#                         If this is 1, use the visualization software to plot
#                         the data
#   --visualization_software VISUALIZATION_SOFTWARE
#                         Use this visualization software to plot the data if
#                         requested
#   --prop_add_var PROP_ADD_VAR
#                         proportion (0-1) of non-flipped regions to have their
#                         methylation calls changed with prob = 0.5, default is
#                         to NO cell to cell variability
#   --bulk_depth BULK_DEPTH
#                         Number of cells that will be used to generate bulk
#                         methylation levels. If zero no bulk data will be
#                         saved.
#   --phylogenetic_generation PHYLOGENETIC_GENERATION
#                         1 or 0. If this is 1, use phylogenetic tree to
#                         generate the clones, if this is 0, the clones are
#                         independent.
#   --prop_cpg_flip PROP_CPG_FLIP
#                         proportion of CpGs to be flipped inside a region
#   --seed SEED           The variability seed. You can set the seed for
#                         reproducibility
#   --verbose VERBOSE     Set to 1 if you want details of the data generation
#   --saveall SAVEALL     Set to 1 if you want the save all the data (with
#                         errors but without missing observations)



