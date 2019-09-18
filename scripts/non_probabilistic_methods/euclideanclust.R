
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

#======================
# EuclideanClust: hierarchical clustering considering Euclidean distances and complete linkage
#======================

if (R == 1){

  if (impute == 1) {
    stop("Imputing option not implemented for this case")
  }

  print("One region, CpG based hiearchical clustering")

  pairwisedist_file <- paste0(outdir, "/pairwisedist.RDa.gz")
  if (file.exists(pairwisedist_file) & use_cache) {
    print("loading previously calculated pairwise dist")
    load(pairwisedist_file)
    print("... done.")
  } else {
    pairwisedist <- dist(input_CpG_data,method="euclidean")
    save(pairwisedist, file = pairwisedist_file, compress = "gzip")
  }

  if(sum(is.na(pairwisedist)) == 0){

    hclust_CpG_crash <- 0

    hcluster <- hclust(pairwisedist,method = "complete")

    # defining some clusters
    mycl <- cutree(hcluster, k=1:Max_K)

    possible_clusters <- cbind(rownames(input_CpG_data),mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id",paste0("EuclideanClust_cpg_num_clusters_",1:Max_K))

    #print(possible_clusters)

    t <- try(NbClust(input_CpG_data, diss = pairwisedist,distance=NULL, min.nc=2, max.nc=Max_K,method = "complete",index = "cindex"))
    if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
      print("can't find a best partition")
      error_ch_index <- 1
    } else {
      error_ch_index <- 0
      hcluster_Nb <- t
      # print(hcluster_Nb)
    }

    write.table(error_ch_index,file=paste0(outdir,"/EuclideanClust_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)

    if(error_ch_index == 0){
      best_cluster <- hcluster_Nb$Best.partition

      possible_clusters <- cbind(possible_clusters,best_cluster)
      colnames(possible_clusters) <- c(colnames(possible_clusters)[1:(dim(possible_clusters)[2]-1)],paste0("best_cluster_",hcluster_Nb$Best.nc[1]))
    }

    ofile <- paste0(outdir,"/EuclideanClust_clusters_CpG_based_maxk_",Max_K,".tsv")
    write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    system(paste("gzip --force", ofile))

    ofile <- paste0(outdir,"/EuclideanClust_cell_order_CpG_based_maxk_",Max_K,".tsv")
    write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste("gzip --force", ofile))
    rm(possible_clusters)

  } else {
    print("some pairs of cells have no CpG with data in common")
    hclust_CpG_crash <- 1
  }
  write.table(hclust_CpG_crash,file=paste0(outdir,"/EuclideanClust_crash.tsv"),row.names=FALSE,col.names=FALSE)
}


if (R > 1){

  print("More than one region, region based hiearchical clustering")

  #======================
  # extracting the mean methylation for each region in each cell
  # this will result in matrix with N cells by R regions - regions are columns and cells the lines
  #======================

  # MA 8 May 2018: saving the mean methylation matrix
  # MA 20 Jan 2019: doing imputation if required
  if (impute == 1) {
    imputed_file <- paste0(outdir,"/EuclideanClust_mean_meth_matrix.Rda.gz")
    if (file.exists(imputed_file) & use_cache) {
      print ("Reading the imputed file")
      load(imputed_file)
      print (" ... done.")
    } else {

      print("Per region, replacing NAs with average values")
      mean_meth_matrix <- impute_means(mean_meth_matrix)
      print(" ... done.")

      # eliminate the empty rows (features)
      mean_meth_matrix <- mean_meth_matrix[ rowSums(mean_meth_matrix)!=0, ]
      save(mean_meth_matrix, file = imputed_file, compress = "gzip")
    }
  }

  pairwisedist_region_file <- paste0(outdir, "/pairwisedist_region.RDa.gz")
  if(file.exists(pairwisedist_region_file) & use_cache){
    print("loading pairwise Euclidean distances")
    load(pairwisedist_region_file)
    print("... done.")
  } else {
    print("Doing hclust on (dis)similarity matrix")
    ### Feb 27th, 2018
    ### doing hclust on the (dis)similarity matrix. Similar to PBAL but based on region mean methylation
    dist_region <- dist(mean_meth_matrix, method="euclidean")
    pairwisedist_region <- dist(dist_region, method="euclidean")
    save(dist_region, pairwisedist_region, file = pairwisedist_region_file, compress = "gzip")
    print("... done.")
  }


  if(sum(is.na(pairwisedist_region)) == 0){

    hclust_region_crash <- 0

    hcluster <- hclust(pairwisedist_region,method = "complete")

    if(!file.exists(paste0(outdir,"/Region_based_EuclideanClust_PLOT.pdf"))){
      pheatmap(as.matrix(dist_region),cluster_rows = TRUE,cluster_cols=TRUE, cellwidth = 8,
               cellheight = 8,fontsize = 8,
               clustering_distance_rows = "euclidean",
               clustering_method = "complete",
               main = paste0("Region-based EuclideanClust"),
               filename = paste0(outdir,"/Region_based_EuclideanClust_PLOT.pdf"))
    }

    # defining some clusters
    mycl <- cutree(hcluster, k=1:Max_K)

    possible_clusters <- cbind(rownames(input_CpG_data),mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id",paste0("EuclideanClust_region_num_clusters_",1:Max_K))

    t <- try(NbClust(as.matrix(dist_region), diss = pairwisedist_region, distance=NULL, min.nc=1, max.nc=Max_K,method = "complete",index = index_type))
    if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
      print("can't use ch or gap index")
      error_ch_index <- 1
    } else {
      error_ch_index <- 0
      hcluster_Nb <- t
      # print(hcluster_Nb)
    }

    hclust_region_bestpartition_crash <- error_ch_index

    write.table(error_ch_index,file=paste0(outdir,"/EuclideanClust_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)

    if(error_ch_index == 0){
      best_cluster <- hcluster_Nb$Best.partition

      possible_clusters <- cbind(possible_clusters,best_cluster)
      colnames(possible_clusters) <- c(colnames(possible_clusters)[1:(dim(possible_clusters)[2]-1)],paste0("best_cluster_",hcluster_Nb$Best.nc[1]))
    }

    ofile <- paste0(outdir,"/EuclideanClust_clusters_region_based_maxk_",Max_K,".tsv")
    write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)

    ofile <- paste0(outdir,"/EuclideanClust_cell_order_region_based_maxk_",Max_K,".tsv")
    write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste("gzip --force", ofile))

    rm(possible_clusters)
  } else {
    print("some pairs of cells have no region with data in common")
    hclust_region_crash <- 1
  }
  write.table(hclust_region_crash,file=paste0(outdir,"/EuclideanClust_crash.tsv"),row.names=FALSE,col.names=FALSE)
}

print("Done Euclidean Clust!")




