
#======================
# libraries
#======================

suppressMessages(library(argparse))

# Renaming the simple methods
# hclust -> EuclideanClust
# PBALclust -> HammingClust
# Pearsonclust -> Pearsonclust
# densitycut -> DensityCut

suppressMessages(library(REpiclomal))

RSCRIPT <- "Rscript"

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--method", default="None", type="character", help="non-probabilistic method to use. densitycut, euclidean, hamming, pearson. all for all methods plus post processing. None for loading and caching data")
parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--methylation_file", type="character", help="Path to methylation data")
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--max_k", type="integer",default=5, help="maximum number of clusters to be considered when cutting the tree")
parser$add_argument("--true_clone_membership_file", default=NULL, type="character",help="Path to true clone membership file")
parser$add_argument("--true_prevalences", default="None", type="character",help="The true prevalence, for example 0.33_0.33_0.34")
parser$add_argument("--evaluate_clustering_software", default="None", type="character",help="Path to the file evaluate_clustering.py")

parser$add_argument("--index", type="character",default="ch",help="Index to be used to choose number of clusters, default = ch")

# MA: 30Jan 2019: adding an optional imputation step
parser$add_argument("--impute", default="0", type="integer",help="If it is 1, impute with the average per region/locus, if it is 0 do nothing.")

parser$add_argument("--use_cache", default="1", type="integer", help="If 1, use cached data if available to save on compute time, if 0, recompute data")

args <- parser$parse_args()

print(args)

method <- args$method
methods <- c('densitycut', 'euclidean', 'hamming', 'pearson', 'all', 'None')
if (!is.element(method, methods)) {
  stop("method must be one of densitycut, euclidean, hamming, pearson, all, None")
}

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)

input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file
true_clusters_file <- args$true_clone_membership_file
eval_soft <- args$evaluate_clustering_software
index_type <- args$index
impute <- args$impute
use_cache <- args$use_cache

## TO TEST
#input_CpG_data_file <- "/Users/camila/Documents/shahlab15/BS-seq/whole_genome_single_cell/synthetic_data_tests/data/data_incomplete.tsv.gz"
#input_regions_file <- "/Users/camila/Documents/shahlab15/BS-seq/whole_genome_single_cell/synthetic_data_tests/data/regions_file.tsv.gz"
#inferred_clusters_file <- "/Users/camila/Documents/shahlab15/BS-seq/whole_genome_single_cell/synthetic_data_tests/data/true_clone_membership.tsv.gz"
#true_clusters_file <- "/Users/camila/Documents/shahlab15/BS-seq/whole_genome_single_cell/synthetic_data_tests/data/true_clone_membership.tsv.gz"
#outdir <- "/Users/camila/Documents/shahlab15/BS-seq/whole_genome_single_cell/synthetic_data_tests/results"
#true_clusters <- read.csv(true_clusters_file ,sep="\t",header=TRUE,check.names=FALSE)


#======================
# loading the data
#======================

# Methylation data
data <- load_data(outdir, input_CpG_data_file, input_regions_file, use_cache)

input_CpG_data <- data$input_CpG_data
input_regions <- data$input_regions
mean_meth_matrix <- data$mean_meth_matrix
rm(data)

R <- dim(input_regions)[1] ## number of regions
M <- dim(input_CpG_data)[2] ## number of loci

Max_K <- min((dim(input_CpG_data)[1]-1),args$max_k)
print(Max_K)

print('number of cells:')
print(dim(input_CpG_data)[1])
print("number of regions:")
print(R)
print("number of loci:")
print(M)

possible_clusters <- NULL

if (method == "densitycut" || method == "all") {
  #################################
  ## DensityCut
  #################################
  # Now call densitycut
  print("Calling DensityCut")

  max_PC <- min(20, R)
  maxit <- 100
  possible_clusters <- densitycut.clust(input_CpG_data, mean_meth_matrix, R, Max_K, max_PC, maxit, impute, use_cache, outdir)
  print("Done DensityCut")
}

if (method == "euclidean" || method == "all") {
  #======================
  # EuclideanClust: hierarchical clustering considering Euclidean distances and complete linkage
  #======================

  # Now call EuclideanClust
  print("Calling EuclideanClust")

  euclidean_clusters <- euclidean.clust(input_CpG_data, mean_meth_matrix, R, Max_K, index_type, impute, use_cache, outdir)

  if (method == "all") {
    possible_clusters <- cbind(possible_clusters, euclidean_clusters[,1:Max_K+1])
    rm(euclidean_clusters)
    hclust_region_crash <- read.table(file=file.path(outdir,"EuclideanClust_crash.tsv"))
    hclust_region_bestpartition_crash <- read.table(file=file.path(outdir,"EuclideanClust_bestpartition_crash.tsv"))
  }
  print("Done Euclidean Clustering")
}

if(method == "hamming" || method == "all") {
  ###################################################
  ### HammingClust: Tony's (PBAL manuscript) clustering approach ##
  ###################################################

  # Now call HammingClust
  print("Calling HammingClust")

  hamming_clusters <- hamming.clust(input_CpG_data, Max_K, index_type, impute, use_cache, outdir)

  if (method == "all") {
    possible_clusters <- cbind(possible_clusters, hamming_clusters[,1:Max_K+1])
    rm(hamming_clusters)
    PBAL_crash <- read.table(file=file.path(outdir,"HammingClust_crash.tsv"))
    PBALclust_bestpartition_crash <- read.table(file=file.path(outdir,"HammingClust_bestpartition_crash.tsv"))
  }
  print("Done Hamming Clustering")
}

if(method == "pearson" || method == "all") {
  ############################################################
  ## PearsonClust: Pearson correlation approach similar to scTrio paper  ###
  ############################################################

  # Now call HammingClust
  print("Calling PearsonClust")

  pearson_clusters <- pearson.clust(input_CpG_data, Max_K, index_type, impute, use_cache, outdir)

  if (method == "all") {
    possible_clusters <- cbind(possible_clusters, pearson_clusters[,1:Max_K+1])
    rm(pearson_clusters)
    Pearson_crash <- read.table(file=file.path(outdir,"PearsonClust_crash.tsv"))
    Pearsonclust_bestpartition_crash <- read.table(file=file.path(outdir,"PearsonClust_bestpartition_crash.tsv"))
  }
  print("Done Pearson Clustering")
}

if (method == "all") {
  #################################
  ## Post Processing
  #################################
  # MA: added another file at the end with all the columns from hclust regions and pbal (except the first 2 columns of pbal cell_id and clusters_1
  # Note: I am unzipping so I zip again after, I should not zip earlier, TODO
  # TODO: skip the best column, this file will be used only for initialization, not for evaluation
  outfile <- gzfile(file.path(outdir, "initial_inputs.tsv.gz"))

  write.table(possible_clusters, file=outfile, col.names=TRUE, sep="\t", quote=FALSE, row.names=FALSE)

  # PYTHON3 <- "/home/mandronescu/.local/centos6/anaconda3/bin/python3"

  # If there is a true clusters file, then run evaluation software
  if (!is.null(true_clusters_file)) {

    if (hclust_region_crash == 0 && hclust_region_bestpartition_crash == 0) {
        print("Calling evaluation software for Hclust (EuclideanClust)")
        command <- paste("python3", eval_soft, "--true_clusters_file", true_clusters_file, "--true_prevalences", args$true_prevalences, "--predicted_clusters_file", paste0(hfile, ".gz"), "--clusters_are_probabilities False --results_file", file.path(outdir, "results_EuclideanClust.txt"))
        print(command)
        system(command)
    }

    if (PBAL_crash == 0 && PBALclust_bestpartition_crash == 0) {
        print("Calling evaluation software for HammingClust")
        command <- paste("python3", eval_soft, "--true_clusters_file", true_clusters_file, "--true_prevalences", args$true_prevalences, "--predicted_clusters_file", paste0(pfile, ".gz"), "--clusters_are_probabilities False --results_file", file.path(outdir, "results_HammingClust.txt"))
        print(command)
        system(command)
    }

    if (Pearson_crash == 0 && Pearsonclust_bestpartition_crash == 0) {
        print("Calling evaluation software for PearsonClust")
        command <- paste("python3", eval_soft, "--true_clusters_file", true_clusters_file, "--true_prevalences", args$true_prevalences, "--predicted_clusters_file", paste0(peafile, ".gz"), "--clusters_are_probabilities False --results_file", file.path(outdir, "results_PearsonClust.txt"))
        print(command)
        system(command)
    }

    if(file.exists(paste0(dfile,".gz"))) {
        print("Calling evaluation software for DensityCut")
        command <- paste("python3", eval_soft, "--true_clusters_file", true_clusters_file, "--true_prevalences", args$true_prevalences, "--predicted_clusters_file", paste0(dfile, ".gz"), "--clusters_are_probabilities False --results_file", file.path(outdir, "results_DensityCut.txt"))
        print(command)
        system(command)
    }
  }

  print("Done!")
}