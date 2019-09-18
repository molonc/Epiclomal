
#======================
# libraries
#======================

suppressMessages(library(argparse))
suppressMessages(library(Rcpp))

# Renaming the simple methods
# hclust -> EuclideanClust
# PBALclust -> HammingClust
# Pearsonclust -> Pearsonclust
# densitycut -> DensityCut


### libraries to find the best number of clusters
#suppressMessages(library(factoextra))
suppressMessages(library(NbClust))
suppressMessages(library(pheatmap))

RSCRIPT <- "Rscript"

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--methylation_file", type="character", help="Path to methylation data")
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--max_k", type="integer",default=5, help="maximum number of clusters to be considered when cutting the tree")
parser$add_argument("--true_clone_membership_file", default=NULL, type="character",help="Path to true clone membership file")
parser$add_argument("--true_prevalences", default="None", type="character",help="The true prevalence, for example 0.33_0.33_0.34")
parser$add_argument("--evaluate_clustering_software", type="character",help="Path to the file evaluate_clustering.py")

parser$add_argument("--index", type="character",default="ch",help="Index to be used to choose number of clusters, default = ch")

# MA: 30Jan 2019: adding an optional imputation step
parser$add_argument("--impute", default="0", type="integer",help="If it is 1, impute with the average per region/locus, if it is 0 do nothing.")

parser$add_argument("--use_cache", default="1", type="integer", help="If 1, use cached data if available to save on compute time, if 0, recompute data")

args <- parser$parse_args()

print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)
input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file
true_clusters_file <- args$true_clone_membership_file
eval_soft <- args$evaluate_clustering_software
impute <- args$impute
use_cache <- args$use_cache

# get directory of current script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(paste(sep="/", script.basename, "helpers.R"))

if(impute == 1){
  sourceCpp(paste(sep="/", script.basename, "cpp_functions", "impute.cpp"))
}


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
data <- load_data(input_CpG_data_file, input_regions_file)

input_CpG_data <- data$input_CpG_data
input_regions <- data$input_regions
mean_meth_matrix <- data$mean_meth_matrix
rm(data)

R <- dim(input_regions)[1] ## number of regions
M <- dim(input_CpG_data)[2] ## number of loci

Max_K <- min((dim(input_CpG_data)[1]-1),args$max_k)
print(Max_K)

index_type <- args$index

#print(R)
#print(M)

#======================
# EuclideanClust: hierarchical clustering considering Euclidean distances and complete linkage
#======================

# Now call EuclideanClust
print("Calling EuclideanClust")

### TO TEST
#stop()

# assume densitycut is in the same directory as this file
# call densitycut.R
euclidean.clust.name <- paste(sep="/", script.basename, "euclideanclust.R")
print(paste("Sourcing",euclidean.clust.name,"from",script.name))

command <- paste(RSCRIPT, euclidean.clust.name,
    "--output_directory", outdir,
    "--max_k", Max_K, "--methylation_file", input_CpG_data_file,
    "--regions_file", input_regions_file,
    "--index", index_type,
    "--impute", impute,
    "--use_cache", use_cache)

print(command)
system(command)

hclust_region_crash <- read.table(file=paste0(outdir,"/EuclideanClust_crash.tsv"))
hclust_region_bestpartition_crash <- read.table(file=paste0(outdir,"/EuclideanClust_bestpartition_crash.tsv"))

###################################################
### HammingClust: Tony's (PBAL manuscript) clustering approach ##
###################################################

# Now call HammingClust
print("Calling HammingClust")

### TO TEST
#stop()

# assume densitycut is in the same directory as this file
# call densitycut.R
hamming.clust.name <- paste(sep="/", script.basename, "hammingclust.R")
print(paste("Sourcing",hamming.clust.name,"from",script.name))

command <- paste(RSCRIPT, hamming.clust.name,
    "--output_directory", outdir,
    "--max_k", Max_K, "--methylation_file", input_CpG_data_file,
    "--regions_file", input_regions_file,
    "--index", index_type,
    "--impute", impute,
    "--use_cache", use_cache)

print(command)
system(command)

PBAL_crash <- read.table(file=paste0(outdir,"/HammingClust_crash.tsv"))
PBALclust_bestpartition_crash <- read.table(file=paste0(outdir,"/HammingClust_bestpartition_crash.tsv"))

############################################################
## PearsonClust: Pearson correlation approach similar to scTrio paper  ###
############################################################

# Now call HammingClust
print("Calling PearsonClust")

### TO TEST
#stop()

# assume densitycut is in the same directory as this file
# call densitycut.R
pearson.clust.name <- paste(sep="/", script.basename, "pearsonclust.R")
print(paste("Sourcing",pearson.clust.name,"from",script.name))

command <- paste(RSCRIPT, pearson.clust.name,
    "--output_directory", outdir,
    "--max_k", Max_K, "--methylation_file", input_CpG_data_file,
    "--regions_file", input_regions_file,
    "--index", index_type,
    "--impute", impute,
    "--use_cache", use_cache)

print(command)
system(command)

Pearson_crash <- read.table(file=paste0(outdir,"/PearsonClust_crash.tsv"))
Pearsonclust_bestpartition_crash <- read.table(file=paste0(outdir,"/PearsonClust_bestpartition_crash.tsv"))


#################################
## DensityCut
#################################
# Now call densitycut
print("Calling DensityCut")

### TO TEST
#stop()

# assume densitycut is in the same directory as this file
# call densitycut.R
densitycut.name <- paste(sep="/", script.basename, "densitycut.R")
print(paste("Sourcing",densitycut.name,"from",script.name))

maxpc <- min(20, R)
command <- paste(RSCRIPT, densitycut.name,
    "--output_directory", outdir,
    "--max_PC", maxpc, "--methylation_file", input_CpG_data_file,
    "--regions_file", input_regions_file,
    "--impute", impute,
    "--use_cache", use_cache)

print(command)
system(command)

#################################
## Post Processing
#################################
# MA: added another file at the end with all the columns from hclust regions and pbal (except the first 2 columns of pbal cell_id and clusters_1
# Note: I am unzipping so I zip again after, I should not zip earlier, TODO
# TODO: skip the best column, this file will be used only for initialization, not for evaluation
if (R == 1) {
  hfile <- paste0(outdir,"/EuclideanClust_clusters_CpG_based_maxk_",Max_K,".tsv")
} else {
  hfile <- paste0(outdir,"/EuclideanClust_clusters_region_based_maxk_",Max_K,".tsv")
}
pfile <- paste0(outdir,"/HammingClust_clusters_CpG_based_maxk_",Max_K,".tsv")
peafile <- paste0(outdir,"/PearsonClust_clusters_CpG_based_maxk_",Max_K,".tsv")
dfile <- paste0(outdir,"/DensityCut_clusters_Region_based_maxPC_",maxpc,".tsv")
#outfile <- paste0(outdir, "/all_hclust_maxk_",Max_K,"_DensityCut_PC", maxpc, ".tsv")
outfile <- paste0(outdir, "/initial_inputs.tsv")
print(paste("Hfile", hfile))
print(paste("Pfile", pfile))
print(paste("Peafile", peafile))
print(paste("Dfile", dfile))
print(paste("Outfile", outfile))

# take the cell_ids from the input methylation file
idtempfile <- paste0(outdir,"/cellid_temp_maxk_",Max_K,".tsv")
htempfile <- ""
ptempfile <- ""
peatempfile <- ""
dtempfile <- ""

#hclust
if(file.exists(hfile)) {
    print(paste0("hclust result exists ", hfile))
    htempfile <- paste0(outdir,"/EuclideanClust_temp_maxk_",Max_K,".tsv")
    system (paste0("cut -f1 ", hfile, " > ", idtempfile))
    system (paste0("cut -f2-", Max_K+1, " ", hfile, " > ", htempfile))
    system (paste0("gzip --force ", hfile))

}

# pbal
if(file.exists(pfile)) {
    print(paste0("HammingClust result exists ", pfile))
    ptempfile <- paste0(outdir,"/HammingClust_temp_maxk_",Max_K,".tsv")
    system (paste0("cut -f1 ", pfile, " > ", idtempfile))
    system (paste0("cut -f2-", Max_K+1, " ", pfile, " > ", ptempfile))
    system (paste0("gzip --force ", pfile))
}

# 12 Apr 2018, adding the Pearsonclust results too
# Pearsonclust
if(file.exists(peafile)) {
    print(paste0("PearsonClust result exists ", peafile))
    peatempfile <- paste0(outdir,"/PearsonClust_temp_maxk_",Max_K,".tsv")
    system (paste0("cut -f1 ", peafile, " > ", idtempfile))
    system (paste0("cut -f2-", Max_K+1, " ", peafile, " > ", peatempfile))
    system (paste0("gzip --force ", peafile))
}

# densitycut
if(file.exists(dfile)) {
    print(paste0("DensityCut result exists ", dfile))
    dtempfile <- paste0(outdir,"/DensityCut_temp_maxpc_",maxpc,".tsv")
    system (paste0("cut -f1 ", dfile, " > ", idtempfile))
    system (paste0("cut -f2 ", dfile, " > ", dtempfile))
    system (paste0("gzip --force ", dfile))
}

 # I should add the name of PBAL or HCLUST - added
# Write the file only when at least one of the files exists
if(file.exists(paste0(hfile,".gz")) || file.exists(paste0(pfile,".gz")) || file.exists(paste0(peafile,".gz")) || file.exists(paste0(dfile,".gz")))  {
    command <- paste0("paste ", idtempfile, " ", htempfile, " ", ptempfile,  " ", peatempfile, " ", dtempfile, " > ", outfile)
    system(command)
    system (paste0("gzip --force ", outfile))
}

if (ptempfile != "")  system (paste("rm", ptempfile))
if (peatempfile != "")  system (paste("rm", peatempfile))
if (htempfile != "")  system (paste("rm", htempfile))
if (dtempfile != "")  system (paste("rm", dtempfile))
if (idtempfile != "")  system (paste("rm", idtempfile))

# PYTHON3 <- "/home/mandronescu/.local/centos6/anaconda3/bin/python3"

# If there is a true clusters file, then run evaluation software
if (!is.null(true_clusters_file)) {

  if (hclust_region_crash == 0 && hclust_region_bestpartition_crash == 0) {
      print("Calling evaluation software for Hclust (EuclideanClust)")
      command <- paste("python3", eval_soft, "--true_clusters_file", true_clusters_file, "--true_prevalences", args$true_prevalences, "--predicted_clusters_file", paste0(hfile, ".gz"), "--clusters_are_probabilities False --results_file", paste0(outdir, "/results_EuclideanClust.txt"))
      print(command)
      system(command)
  }

  if (PBAL_crash == 0 && PBALclust_bestpartition_crash == 0) {
      print("Calling evaluation software for HammingClust")
      command <- paste("python3", eval_soft, "--true_clusters_file", true_clusters_file, "--true_prevalences", args$true_prevalences, "--predicted_clusters_file", paste0(pfile, ".gz"), "--clusters_are_probabilities False --results_file", paste0(outdir, "/results_HammingClust.txt"))
      print(command)
      system(command)
  }

  if (Pearson_crash == 0 && Pearsonclust_bestpartition_crash == 0) {
      print("Calling evaluation software for PearsonClust")
      command <- paste("python3", eval_soft, "--true_clusters_file", true_clusters_file, "--true_prevalences", args$true_prevalences, "--predicted_clusters_file", paste0(peafile, ".gz"), "--clusters_are_probabilities False --results_file", paste0(outdir, "/results_PearsonClust.txt"))
      print(command)
      system(command)
  }

  if(file.exists(paste0(dfile,".gz"))) {
      print("Calling evaluation software for DensityCut")
      command <- paste("python3", eval_soft, "--true_clusters_file", true_clusters_file, "--true_prevalences", args$true_prevalences, "--predicted_clusters_file", paste0(dfile, ".gz"), "--clusters_are_probabilities False --results_file", paste0(outdir, "/results_DensityCut.txt"))
      print(command)
      system(command)
  }
}

print("Done!")
