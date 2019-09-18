
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

############################################################
## PearsonClust: Pearson correlation approach similar to scTrio paper  ###
############################################################

dist_Pearson_file <- paste0(outdir, "/dist_Pearson.RDa.gz")
if (file.exists(dist_Pearson_file) & use_cache) {
  print("Loading pearson distances from file")
  load(dist_Pearson_file)
  print("... done.")
} else {
  print("computing Pearson distances")
  dist_Pearson <- cor(x=t(input_CpG_data),method="pearson", use ="pairwise.complete.obs")
  ### from PBAL manuscript: unsupervised learning was done by calculating a Euclidean distance from each cell’s dissimilarity vector and clustered using Ward’s linkage method.
  print("scTrio's approach - CpG based clustering")
  diss_matrix_P <- 1 - cor(x=dist_Pearson,method="pearson", use ="pairwise.complete.obs")
  print("done computing, saving to file")
  save(dist_Pearson, diss_matrix_P, file = dist_Pearson_file, compress = "gzip")
}

if(sum(is.na(diss_matrix_P)) == 0){

  Pearson_crash <- 0

  hcluster_P <- hclust(as.dist(diss_matrix_P),method = "ward.D2") ### BECAUSE IT IS CORRELATION diss_matrix_P IS ACTUALLY A SIMILARITY MATRIX, SO HAVE TO DO 1 - diss_matrix_P

  if(!file.exists(paste0(outdir,"/PearsonClust_PLOT.pdf"))){
    pheatmap(dist_Pearson ,cluster_rows = TRUE,cluster_cols=TRUE, cellwidth = 8,
             cellheight = 8,fontsize = 8,
             #clustering_distance_rows = "correlation",
             #clustering_distance_cols = "correlation",
             clustering_distance_cols = as.dist(diss_matrix_P),
             clustering_distance_rows = as.dist(diss_matrix_P),
             clustering_method = "ward.D2",
             main = paste0("Pearson corr. approach"),
             filename = paste0(outdir,"/PearsonClust_PLOT.pdf"))
  }

  # defining some clusters
  mycl_P <- cutree(hcluster_P, k=1:Max_K)

  possible_clusters_P <- cbind(rownames(input_CpG_data),mycl_P)
  possible_clusters_P <- as.data.frame(possible_clusters_P)
  colnames(possible_clusters_P) <- c("cell_id",paste0("PearsonClust_num_clusters_",1:Max_K))

  t <- try(NbClust(dist_Pearson, diss = as.dist(diss_matrix_P),distance=NULL, min.nc=1, max.nc=Max_K,method = "ward.D2",index = index_type))
  if("try-error" %in% class(t)) {
    print("can't find best partition")
    error_ch_index <- 1
  } else {
      error_ch_index <- 0
      hcluster_Nb_P <- NbClust(dist_Pearson , diss = as.dist(diss_matrix_P),distance=NULL, min.nc=1, max.nc=Max_K,method = "ward.D2",index = index_type)
      # print(hcluster_Nb_P)
    }

  Pearsonclust_bestpartition_crash <- error_ch_index

  write.table(error_ch_index,file=paste0(outdir,"/PearsonClust_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)

  if(error_ch_index == 0){
    best_cluster <- hcluster_Nb_P$Best.partition

    possible_clusters_P <- cbind(possible_clusters_P,best_cluster)
    colnames(possible_clusters_P) <- c(colnames(possible_clusters_P)[1:(dim(possible_clusters_P)[2]-1)],paste0("best_cluster_",hcluster_Nb_P$Best.nc[1]))
  }

  ofile <- paste0(outdir,"/PearsonClust_clusters_CpG_based_maxk_",Max_K,".tsv")
  write.table(possible_clusters_P, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)

  ofile <- paste0(outdir,"/PearsonClust_cell_order_CpG_based_maxk_",Max_K,".tsv")
  write.table(hcluster_P$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
  system(paste0("gzip --force ", ofile))

  rm(possible_clusters_P)
} else {
  Pearson_crash <- 1
}
write.table(Pearson_crash,file=paste0(outdir,"/PearsonClust_crash.tsv"),row.names=FALSE,col.names=FALSE)

print("Done Pearson Clust!")




