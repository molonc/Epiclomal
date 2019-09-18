
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

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--methylation_file", type="character", help="Path to methylation data")
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--max_k", type="integer",default=5, help="maximum number of clusters to be considered when cutting the tree")

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

###################################################
### HammingClust: Tony's (PBAL manuscript) clustering approach ##
###################################################

### from PBAL manuscript: unsupervised learning was done by calculating a Euclidean distance from each cell’s dissimilarity vector and clustered using Ward’s linkage method.

print("Tony's approach - CpG based clustering (HammingClust)")

#input_CpG_data <- input_CpG_data[1:100,1:5]
if (impute == 1) {
  imputed_file <- paste0(outdir,"/CpG_based_imputed.RDa.gz")
  if (file.exists(imputed_file) & use_cache) {
    print ("Reading the imputed file")
    load(imputed_file)
    print (" ... done.")
  } else {
    # replace with average values, for each col
    print("Per region, replacing NAs with median values")
    input_CpG_data <- impute_medians(input_CpG_data)
    print(" ... done.")

    # eliminate the empty rows (features)
    input_CpG_data <- input_CpG_data[ rowSums(input_CpG_data)!=0, ]
    save(input_CpG_data, file = imputed_file, compress = "gzip")
  }
}

dist_PBAL_file <- paste0(outdir, "/dist_PBAL.RDa.gz")
if (file.exists(dist_PBAL_file) & use_cache) {
  print("Loading PBAL distance matrix from file")
  load(dist_PBAL_file)
  print("... done.")
} else {
  # assume dist_PBAL is in the same directory as this file
  # trying to figure out the path of this file so I can call dist_PBAL.cpp
  print(paste("getting dist_pbal from file", paste(sep="/", script.basename, "dist_PBAL.cpp")))
  sourceCpp(paste(sep="/", script.basename, "cpp_functions", "dist_PBAL.cpp"))
  print('Computing PBAL distance matrix')
  dist_PBAL <- as.dist(dist_PBAL(d = input_CpG_data))
  print('Computing pairwise PBAL distance matrix')
  diss_matrix_T <- dist(dist_PBAL,method="euclidean")
  print("Done Computing, saving to file")
  save(dist_PBAL, diss_matrix_T, file = dist_PBAL_file, compress = "gzip")
}

# dist_PBAL <- dist.PBAL(d=input_CpG_data)

# diss_matrix_T <- dist(dist_PBAL,method="euclidean")

if(sum(is.na(diss_matrix_T)) == 0){

  PBAL_crash <- 0

  hcluster_T <- hclust(diss_matrix_T,method = "ward.D2")

  if(!file.exists(paste0(outdir,"/HammingClust_PLOT.pdf"))){
    pheatmap(dist_PBAL,cluster_rows = TRUE,cluster_cols=TRUE, cellwidth = 8,
             cellheight = 8,fontsize = 8,
             clustering_distance_rows = "euclidean",
             clustering_method = "ward.D2",
             main = paste0("HammingClust"),
             filename = paste0(outdir,"/HammingClust_PLOT.pdf"))
  }

  ## defining some clusters
  mycl_T <- cutree(hcluster_T, k=1:Max_K)

  possible_clusters_T <- cbind(rownames(input_CpG_data),mycl_T)
  possible_clusters_T <- as.data.frame(possible_clusters_T)
  colnames(possible_clusters_T) <- c("cell_id",paste0("HammingClust_num_clusters_",1:Max_K))

  t <- try(NbClust(dist_PBAL, diss = diss_matrix_T,distance=NULL, min.nc=1, max.nc=Max_K,method = "ward.D2",index = index_type)) ### changing to cindex as cindex also works for CpG based clustering
  if("try-error" %in% class(t)) {
    print("can't find best partition")
    error_ch_index <- 1
  } else {
    error_ch_index <- 0
    hcluster_Nb_T <- t
    # print(hcluster_Nb_T)
  }

  PBALclust_bestpartition_crash <- error_ch_index

  write.table(error_ch_index,file=paste0(outdir,"/HammingClust_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)

  if(error_ch_index == 0){

    best_cluster <- hcluster_Nb_T$Best.partition

    possible_clusters_T <- cbind(possible_clusters_T,best_cluster)
    colnames(possible_clusters_T) <- c(colnames(possible_clusters_T)[1:(dim(possible_clusters_T)[2]-1)],paste0("best_cluster_",hcluster_Nb_T$Best.nc[1]))
  }

  ofile <- paste0(outdir,"/HammingClust_clusters_CpG_based_maxk_",Max_K,".tsv")
  write.table(possible_clusters_T, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)

  ofile <- paste0(outdir,"/HammingClust_cell_order_CpG_based_maxk_",Max_K,".tsv")
  write.table(hcluster_T$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
  system(paste("gzip --force", ofile))

  rm(possible_clusters_T)
} else {
  PBAL_crash <- 1
}
write.table(PBAL_crash,file=paste0(outdir,"/HammingClust_crash.tsv"),row.names=FALSE,col.names=FALSE)

print("Done Hamming Clust!")




