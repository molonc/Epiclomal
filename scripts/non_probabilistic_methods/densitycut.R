
#======================
# libraries
#======================

### for testing I used tehse paths for the libraries and R in /gsc/software/linux-x86_64-centos6/R-3.3.2/bin/R
# library(pcaMethods,lib.loc="/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.2")
# library(densitycut,lib.loc = "/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.3/" )

suppressMessages(library(argparse))
suppressMessages(library(REpiclomal))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--methylation_file", type="character", help="Path to methylation data")
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--max_PC", type="integer",default=20, help="maximum number of PC components")
parser$add_argument("--maxit", type="integer",default=100, help="maximum number iterations")
# MA: 30Jan 2019: adding an optional imputation step
parser$add_argument("--impute", default="0", type="integer",help="If it is 1, impute with the average per region/locus, if it is 0 do nothing.")
parser$add_argument("--use_cache", default="1", type="integer", help="If 1, use cached data if available to save on compute time, if 0, recompute data")

args <- parser$parse_args()

print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)
input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file
max_PC <- args$max_PC
maxit <- args$maxit
impute <- args$impute
use_cache <- args$use_cache

# input_CpG_data_file <- "data_incomplete.tsv.gz"
# input_regions_file <- "regions_file.tsv.gz"

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# #outdir <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data"
#
#  input_CpG_data_file <- "/genesis/shahlab/csouza/BS-seq/whole_genome_single_cell/synthetic_data_tests/data_old_way/data/data_incomplete.tsv.gz"
#  input_regions_file <- "/genesis/shahlab/csouza/BS-seq/whole_genome_single_cell/synthetic_data_tests/data_old_way/data/regions_file.tsv.gz"
#  inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform_seed_2/data/true_clone_membership.tsv"
#  true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform_seed_2/data/true_clone_membership.tsv"

#  input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
#  input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
#  inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
#  true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci5000_clones3_cells100_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.85_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci5000_clones3_cells100_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.85_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci5000_clones3_cells100_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.85_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci5000_clones3_cells100_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.85_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"

#input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells50_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.95_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
#input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells50_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.95_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
#inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells50_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.95_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
#true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells50_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.95_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"


#======================
# loading the data
#======================

data <- load_data(outdir, input_CpG_data_file, input_regions_file, use_cache)

input_CpG_data <- data$input_CpG_data
input_regions <- data$input_regions
mean_meth_matrix <- data$mean_meth_matrix
rm(data)

R <- dim(input_regions)[1] ## number of regions
M <- dim(input_CpG_data)[2] ## number of loci

print("number of regions:")
print(R)
print("number of loci:")
print(M)

# ########################################
# ### Some tests with bivariate normal ###
# ########################################
#
# #Sigma <- matrix(c(10,3,3,2),nrow=2,ncol=2)
# Sigma <- matrix(c(1,0,0,1),nrow=2,ncol=2)
# Sigma_2 <- matrix(c(1,0,0,1),nrow=2,ncol=2)
#
# test <- mvrnorm(n = 100, rep(0, 2), Sigma)
# test2 <- mvrnorm(n = 100, rep(10, 2), Sigma_2)
#
# data <- rbind(test,test2)
# data2 <- t(apply(data,1,function(x){(x-colMeans(data))}))
#
# pc <- pca(data2,method="nipals",nPcs=2,center=FALSE,scale=NULL)
#
# tmp <- loadings(pc)
#
# tmp2 <- scores(pc)
#
# tmp[,1]%*%data2[1,] == tmp2[1,1]
#
# slplot(pc)

#======================
# hiearchical clustering considering Euclidean distances and complete linkage
#======================
possible_clusters <- densitycut.clust(input_CpG_data, mean_meth_matrix, R, max_PC, maxit, impute, use_cache, outdir)

print("DONE!")




