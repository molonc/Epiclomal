
#======================
# libraries
#======================

suppressMessages(library(argparse))

### libraries to find the best number of clusters
#suppressMessages(library(factoextra))
suppressMessages(library(NbClust))
suppressMessages(library(pheatmap))

RSCRIPT <- "/gsc/software/linux-x86_64-centos6/R-3.3.2/bin/Rscript"

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
parser$add_argument("--true_prevalences", default=NULL, type="character",help="The true prevalence, for example 0.33_0.33_0.34")
parser$add_argument("--evaluate_clustering_software", type="character",help="Path to the file evaluate_clustering.py")

parser$add_argument("--index", type="character",default="ch",help="Index to be used to choose number of clusters, default = ch")

args <- parser$parse_args() 

print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)
input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file
true_clusters_file <- args$true_clone_membership_file
eval_soft <- args$evaluate_clustering_software

#input_CpG_data_file <- "/Users/cdesouza/Documents/shahlab15/BS-seq/whole_genome_single_cell/synthetic_data_tests/data/data_incomplete.tsv.gz"
#input_regions_file <- "/Users/cdesouza/Documents/shahlab15/BS-seq/whole_genome_single_cell/synthetic_data_tests/data/regions_file.tsv.gz"
#inferred_clusters_file <- "/Users/cdesouza/Documents/shahlab15/BS-seq/whole_genome_single_cell/synthetic_data_tests/data/true_clone_membership.tsv.gz"
#true_clusters_file <- "/Users/cdesouza/Documents/shahlab15/BS-seq/whole_genome_single_cell/synthetic_data_tests/data/true_clone_membership.tsv.gz"

#======================
# auxiliar functions
#======================

extract_mean_meth_per_cell <- function(cell_data,region_coord){
  mean_meth <- apply(region_coord,1,function(x){mean(cell_data[x[1]:x[2]],na.rm=TRUE)})  
  mean_meth[is.na(mean_meth)] <- NA
  return(mean_meth)
}  

#======================
# loading the data
#======================

# Methylation data
tmp <- read.csv(input_CpG_data_file,sep="\t",header=TRUE,check.names=FALSE)
input_CpG_data <- as.matrix(tmp[,-1])
rownames(input_CpG_data) <- tmp$cell_id
rm(tmp)

# Region coordinates
tmp <- read.csv(input_regions_file,sep="\t",header=TRUE,check.names=FALSE)
input_regions <- as.matrix(tmp[,-1]) + 1 ## adding 1 to match R indexing - previously coordinates were for python starting on zero
colnames(input_regions) <- c("start","end") ## input_regions gives already the columns in input_CpG_data that correspond to which regions considered in the construction of input_CpG_data
rownames(input_regions) <- tmp$region_id
rm(tmp)

R <- dim(input_regions)[1] ## number of regions
M <- dim(input_CpG_data)[2] ## number of loci

Max_K <- min((dim(input_CpG_data)[1]-1),args$max_k)
print(Max_K)

index_type <- args$index

#print(R)
#print(M)

#======================
# hiearchical clustering considering Euclidean distances and complete linkage 
#======================

if (R == 1){
  
  print("One region, CpG based hiearchical clustering")
  
  pairwisedist <- dist(input_CpG_data,method="euclidean")
  
  if(sum(is.na(pairwisedist)==TRUE) == 0){
    
    hclust_CpG_crash <- 0
    write.table(hclust_CpG_crash,file=paste0(outdir,"/hclust_CpG_crash.tsv"),row.names=FALSE,col.names=FALSE)
    
    hcluster <- hclust(pairwisedist,method = "complete")
    
    # defining some clusters
    mycl <- cutree(hcluster, k=1:Max_K)
    
    possible_clusters <- cbind(rownames(input_CpG_data),mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id",paste0("hclust_cpg_num_clusters_",1:Max_K))
    
    #print(possible_clusters)
    
    t <- try(NbClust(input_CpG_data, diss = pairwisedist,distance=NULL, min.nc=2, max.nc=Max_K,method = "complete",index = "cindex"))
    if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
      print("can't find a best partition")
      error_ch_index <- 1 }
    else {
      error_ch_index <- 0
      hcluster_Nb <- NbClust(input_CpG_data, diss = pairwisedist,distance=NULL, min.nc=2, max.nc=Max_K,method = "complete",index = "cindex")
      print(hcluster_Nb)}
    if(error_ch_index == 1){
      write.table(error_ch_index,file=paste0(outdir,"/hclust_CpGbased_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)
      
      ofile <- paste0(outdir,"/hclust_clusters_CpG_based_maxk_",Max_K,".tsv")
      write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
      system(paste0("gzip --force ", ofile))
      
      ofile <- paste0(outdir,"/hclust_cell_order_CpG_based_maxk_",Max_K,".tsv")
      write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
      system(paste0("gzip --force ", ofile))
    }
    
    if(error_ch_index == 0){
      
      write.table(error_ch_index,file=paste0(outdir,"/hclust_CpGbased_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)
      
      best_cluster <- hcluster_Nb$Best.partition
      
      possible_clusters <- cbind(possible_clusters,best_cluster)
      colnames(possible_clusters) <- c(colnames(possible_clusters)[1:(dim(possible_clusters)[2]-1)],paste0("best_cluster_",hcluster_Nb$Best.nc[1]))
      
      ofile <- paste0(outdir,"/hclust_clusters_CpG_based_maxk_",Max_K,".tsv") 
      write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
      system(paste0("gzip --force ", ofile))    
      
      ofile <- paste0(outdir,"/hclust_cell_order_CpG_based_maxk_",Max_K,".tsv") 
      write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
      system(paste0("gzip --force ", ofile))    
      
      
    }
    
    rm(possible_clusters)
    
  }
  
  if(sum(is.na(pairwisedist)==TRUE) > 0){
    print("some pairs of cells have no CpG with data in common")
    
    hclust_CpG_crash <- 1
    write.table(hclust_CpG_crash,file=paste0(outdir,"/hclust_CpG_crash.tsv"),row.names=FALSE,col.names=FALSE)
    # system(paste0("gzip --force ", paste0(outdir,"/hclust_crash.tsv"))) 
  }
  
}


if (R > 1){
  
  print("More than one region, region based hiearchical clustering")
  
  #======================
  # extracting the mean methylation for each region in each cell 
  # this will result in matrix with N cells by R regions - regions are columns and cells the lines 
  #======================
  
  mean_meth_matrix <- t(apply(input_CpG_data,1,extract_mean_meth_per_cell,region_coord=input_regions))
  
  # pairwisedist_region <- dist(mean_meth_matrix ,method="euclidean")
  
  ### Feb 27th, 2018
  ### doing hclust on the (dis)similarity matrix. Similar to PBAL but based on region mean methylation
  pairwisedist_region <- dist(dist(mean_meth_matrix ,method="euclidean"),method="euclidean")
  
  if(sum(is.na(pairwisedist_region)==TRUE) == 0){
    
    hclust_region_crash <- 0
    write.table(hclust_region_crash,file=paste0(outdir,"/hclust_region_crash.tsv"),row.names=FALSE,col.names=FALSE)
    
    hcluster <- hclust(pairwisedist_region,method = "complete")
    
    pheatmap(as.matrix(dist(mean_meth_matrix ,method="euclidean")),cluster_rows = TRUE,cluster_cols=TRUE, cellwidth = 8,
             cellheight = 8,fontsize = 8, 
             clustering_distance_rows = "euclidean",
             clustering_method = "complete",
             main = paste0("Region-based hclust"),
             filename = paste0(outdir,"/Region_based_hclust_PLOT.pdf"))
    
    # defining some clusters
    mycl <- cutree(hcluster, k=1:Max_K)
    
    possible_clusters <- cbind(rownames(input_CpG_data),mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id",paste0("hclust_region_num_clusters_",1:Max_K))

    t <- try(NbClust(as.matrix(dist(mean_meth_matrix, method="euclidean")), diss = pairwisedist_region, distance=NULL, min.nc=1, max.nc=Max_K,method = "complete",index = index_type)) 
    if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
      print("can't use ch or gap index") 
      error_ch_index <- 1 } else { 
      error_ch_index <- 0
      hcluster_Nb <- NbClust(as.matrix(dist(mean_meth_matrix, method="euclidean")), diss = pairwisedist_region, distance=NULL, min.nc=1, max.nc=Max_K,method = "complete",index = index_type)
      print(hcluster_Nb)
    }
    
    hclust_region_bestpartition_crash <- error_ch_index
    
    
    if(error_ch_index == 1){
      write.table(error_ch_index,file=paste0(outdir,"/hclust_region_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)
      
      ofile <- paste0(outdir,"/hclust_clusters_region_based_maxk_",Max_K,".tsv")
      write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
      system(paste0("gzip --force ", ofile))
      
      ofile <- paste0(outdir,"/hclust_cell_order_region_based_maxk_",Max_K,".tsv")
      write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
      system(paste0("gzip --force ", ofile))
    }
    
    if(error_ch_index == 0){
      
      write.table(error_ch_index,file=paste0(outdir,"/hclust_region_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)
      
      best_cluster <- hcluster_Nb$Best.partition
      
      possible_clusters <- cbind(possible_clusters,best_cluster)
      colnames(possible_clusters) <- c(colnames(possible_clusters)[1:(dim(possible_clusters)[2]-1)],paste0("best_cluster_",hcluster_Nb$Best.nc[1]))
      
      ofile <- paste0(outdir,"/hclust_clusters_region_based_maxk_",Max_K,".tsv")
      write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
      system(paste0("gzip --force ", ofile))
      
      ofile <- paste0(outdir,"/hclust_cell_order_region_based_maxk_",Max_K,".tsv")
      write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
      system(paste0("gzip --force ", ofile))
      
    }
    
    rm(possible_clusters)
    
  }
  
  if(sum(is.na(pairwisedist_region)==TRUE) > 0){
    print("some pairs of cells have no region with data in common")
    hclust_region_crash <- 1
    write.table(hclust_region_crash,file=paste0(outdir,"/hclust_region_crash.tsv"),row.names=FALSE,col.names=FALSE)
  }
  
}


###################################################
### Tony's (PBAL manuscript) clustering approach ##
###################################################

dist.pair <- function(v1,v2){
  na.idx <- is.na(v1) | is.na(v2) 
  v1a  <- v1[!na.idx]
  v2a  <- v2[!na.idx]
  l.na <- (sum(!na.idx)) ## = length(v1a) = length(v2a), the number of entries with data on both vectors
  
  d <- (sum(abs(v1a - v2a)) / l.na) ### this should be what dist() does! 
  return(d)
  
}

dist.PBAL <- function(d){ ### d is matrix where the rows correspond to cells and columns to CpGs
  dist.data <- matrix(NA,nrow=nrow(d),ncol=nrow(d))
  rownames(dist.data) <- rownames(d)
  colnames(dist.data) <- rownames(d)
  for (i in 1:nrow(d)){
    for(j in 1:nrow(d)){
      dist.data[i,j] <- dist.pair(v1=d[i,],v2=d[j,])
    }
  }
  return(dist.data)
}

### from PBAL manuscript: unsupervised learning was done by calculating a Euclidean distance from each cell’s dissimilarity vector and clustered using Ward’s linkage method.

print("Tony's approach - CpG based clustering")

dist_PBAL <- dist.PBAL(d=input_CpG_data)

diss_matrix_T <- dist(dist_PBAL,method="euclidean")

if(sum(is.na(diss_matrix_T)) == 0){
  
  PBAL_crash <- 0
  write.table(PBAL_crash,file=paste0(outdir,"/PBAL_crash.tsv"),row.names=FALSE,col.names=FALSE)
  
  hcluster_T <- hclust(diss_matrix_T,method = "ward.D2")
  
  pheatmap(dist_PBAL,cluster_rows = TRUE,cluster_cols=TRUE, cellwidth = 8,
           cellheight = 8,fontsize = 8, 
           clustering_distance_rows = "euclidean",
           clustering_method = "ward.D2",
           main = paste0("PBAL approach"),
           filename = paste0(outdir,"/PBAL_hclust_PLOT.pdf"))
  
  ## defining some clusters
  mycl_T <- cutree(hcluster_T, k=1:Max_K)
  
  possible_clusters_T <- cbind(rownames(input_CpG_data),mycl_T)
  possible_clusters_T <- as.data.frame(possible_clusters_T)
  colnames(possible_clusters_T) <- c("cell_id",paste0("pbal_num_clusters_",1:Max_K))
  
  t <- try(NbClust(dist_PBAL, diss = diss_matrix_T,distance=NULL, min.nc=1, max.nc=Max_K,method = "ward.D2",index = index_type)) ### changing to cindex as cindex also works for CpG based clustering
  if("try-error" %in% class(t)) { 
    print("can't find best partition") 
    error_ch_index <- 1 }
  else { 
    error_ch_index <- 0
    hcluster_Nb_T <- NbClust(dist_PBAL, diss = diss_matrix_T,distance=NULL, min.nc=1, max.nc=Max_K,method = "ward.D2",index = index_type)
    print(hcluster_Nb_T)
  }
  
  PBALclust_bestpartition_crash <- error_ch_index
  
  if(error_ch_index == 1){
    write.table(error_ch_index,file=paste0(outdir,"/PBALclust_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)
    
    ofile <- paste0(outdir,"/PBALclust_clusters_CpG_based_maxk_",Max_K,".tsv") 
    write.table(possible_clusters_T, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    system(paste0("gzip --force ", ofile))
    
    ofile <- paste0(outdir,"/PBALclust_cell_order_CpG_based_maxk_",Max_K,".tsv")
    write.table(hcluster_T$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste0("gzip --force ", ofile))
  }
  
  if(error_ch_index == 0){
    
    write.table(error_ch_index,file=paste0(outdir,"/PBALclust_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)
    
    best_cluster <- hcluster_Nb_T$Best.partition
    
    possible_clusters_T <- cbind(possible_clusters_T,best_cluster)
    colnames(possible_clusters_T) <- c(colnames(possible_clusters_T)[1:(dim(possible_clusters_T)[2]-1)],paste0("best_cluster_",hcluster_Nb_T$Best.nc[1]))
    
    ofile <- paste0(outdir,"/PBALclust_clusters_CpG_based_maxk_",Max_K,".tsv") 
    write.table(possible_clusters_T, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    system(paste0("gzip --force ", ofile))
    
    ofile <- paste0(outdir,"/PBALclust_cell_order_CpG_based_maxk_",Max_K,".tsv")
    write.table(hcluster_T$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste0("gzip --force ", ofile))
    
  }
}

if(sum(is.na(diss_matrix_T)) > 0){
  
  PBAL_crash <- 1
  write.table(PBAL_crash,file=paste0(outdir,"/PBAL_crash.tsv"),row.names=FALSE,col.names=FALSE)}


############################################################
## Pearson correlation approach similar to scTrio paper  ###
############################################################

dist.corr <- function(v1,v2){
  na.idx <- is.na(v1) | is.na(v2) 
  v1a  <- v1[!na.idx]
  v2a  <- v2[!na.idx]
  d <- cor(x=v1a,y=v2a)
  return(d)
  
}

dist_Pearson <- cor(x=t(input_CpG_data),method="pearson", use ="pairwise.complete.obs")

# d <- dist.corr(v1=input_CpG_data[1,],v2=input_CpG_data[2,]) ### my own way same as above 

### from PBAL manuscript: unsupervised learning was done by calculating a Euclidean distance from each cell’s dissimilarity vector and clustered using Ward’s linkage method.

print("scTrio's approach - CpG based clustering")

diss_matrix_P <- cor(x=dist_Pearson,method="pearson", use ="pairwise.complete.obs")

if(sum(is.na(diss_matrix_P)) == 0){

  Pearson_crash <- 0
  write.table(Pearson_crash,file=paste0(outdir,"/Pearson_crash.tsv"),row.names=FALSE,col.names=FALSE)

  hcluster_P <- hclust(as.dist(1-diss_matrix_P),method = "ward.D2") ### BECAUSE IT IS CORRELATION diss_matrix_P IS ACTUALLY A SIMILARITY MATRIX, SO HAVE TO DO 1 - diss_matrix_P

  pheatmap(dist_Pearson ,cluster_rows = TRUE,cluster_cols=TRUE, cellwidth = 8,
           cellheight = 8,fontsize = 8,
           clustering_distance_rows = "correlation",
           #clustering_distance_rows = as.dist(1-diss_matrix_P),
           clustering_method = "ward.D2",
           main = paste0("Pearson corr. approach"),
           filename = paste0(outdir,"/Pearson_hclust_PLOT.pdf"))

  # defining some clusters
  mycl_P <- cutree(hcluster_P, k=1:Max_K)

  possible_clusters_P <- cbind(rownames(input_CpG_data),mycl_P)
  possible_clusters_P <- as.data.frame(possible_clusters_P)
  colnames(possible_clusters_P) <- c("cell_id",paste0("pearson_num_clusters_",1:Max_K))

  t <- try(NbClust(dist_Pearson , diss = as.dist(1-diss_matrix_P),distance=NULL, min.nc=1, max.nc=Max_K,method = "ward.D2",index = index_type)) 
  if("try-error" %in% class(t)) { 
    print("can't find best partition")
    error_ch_index <- 1 } else {
    error_ch_index <- 0
    hcluster_Nb_P <- NbClust(dist_Pearson , diss = as.dist(1-diss_matrix_P),distance=NULL, min.nc=1, max.nc=Max_K,method = "ward.D2",index = index_type)
    print(hcluster_Nb_P)
  }

 Pearsonclust_bestpartition_crash <- error_ch_index

  if(error_ch_index == 1){
    write.table(error_ch_index,file=paste0(outdir,"/Pearsonclust_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)

    ofile <- paste0(outdir,"/Pearsonclust_clusters_CpG_based_maxk_",Max_K,".tsv")
    write.table(possible_clusters_P, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    system(paste0("gzip --force ", ofile))

    ofile <- paste0(outdir,"/Pearsonclust_cell_order_CpG_based_maxk_",Max_K,".tsv")
    write.table(hcluster_P$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste0("gzip --force ", ofile))
  }

  if(error_ch_index == 0){

    write.table(error_ch_index,file=paste0(outdir,"/Pearsonclust_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)

    best_cluster <- hcluster_Nb_P$Best.partition

    possible_clusters_P <- cbind(possible_clusters_P,best_cluster)
    colnames(possible_clusters_P) <- c(colnames(possible_clusters_P)[1:(dim(possible_clusters_P)[2]-1)],paste0("best_cluster_",hcluster_Nb_P$Best.nc[1]))

    ofile <- paste0(outdir,"/Pearsonclust_clusters_CpG_based_maxk_",Max_K,".tsv")
    write.table(possible_clusters_P, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    system(paste0("gzip --force ", ofile))

    ofile <- paste0(outdir,"/Pearsonclust_cell_order_CpG_based_maxk_",Max_K,".tsv")
    write.table(hcluster_P$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste0("gzip --force ", ofile))

  }
}

if(sum(is.na(diss_matrix_P)) > 0){

  Pearson_crash <- 1
  write.table(Pearson_crash,file=paste0(outdir,"/Pearson_crash.tsv"),row.names=FALSE,col.names=FALSE)}


####################
# Now call densitycut
print("Calling densitycut")

# assume densitycut is in the same directory as this file
# trying to figure out the path of this file so I can call densitycut.R
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
densitycut.name <- paste(sep="/", script.basename, "densitycut.R")
# print(paste("Sourcing",densitycut.name,"from",script.name))

maxpc <- 20
command <- paste0(RSCRIPT, " ", densitycut.name, 
    " --output_directory ", outdir, 
    " --max_PC ", maxpc, " --methylation_file ", input_CpG_data_file,
    " --regions_file ", input_regions_file)

print(command)
system(command)


# MA: added another file at the end with all the columns from hclust regions and pbal (except the first 2 columns of pbal cell_id and clusters_1
# Note: I am unzipping so I zip again after, I should not zip earlier, TODO
# TODO: skip the best column, this file will be used only for initialization, not for evaluation

hfile <- paste0(outdir,"/hclust_clusters_region_based_maxk_",Max_K,".tsv")
pfile <- paste0(outdir,"/PBALclust_clusters_CpG_based_maxk_",Max_K,".tsv")
peafile <- paste0(outdir,"Pearsonclust_clusters_CpG_based_maxk_",Max_K,".tsv")
dfile <- paste0(outdir,"/densityCut_clusters_Region_based_maxPC_",maxpc,".tsv")
outfile <- paste0(outdir, "/hclust_region_PBAL_CpG_clusters_maxk_",Max_K,"_densitycut_PC", maxpc, ".tsv")
print (paste0("Hfile ", hfile))
print (paste0("Pfile ", pfile))
print (paste0("Dfile ", dfile))
print (paste0("Outfile ", outfile))

# take the cell_ids from the input methylation file
idtempfile <- paste0(outdir,"/cellid_temp_maxk_",Max_K,".tsv")
htempfile <- ""
ptempfile <- ""
dtempfile <- ""

#hclust
if(file.exists(paste0(hfile,".gz"))) {
    print(paste0("hclust result exists ", hfile, ".gz"))
    system(paste0 ("gunzip ", hfile,".gz"))
    htempfile <- paste0(outdir,"/hclust_temp_maxk_",Max_K,".tsv")
    system (paste0 ("cut -f1 ", hfile, " > ", idtempfile))
    system (paste0 ("cut -f2-11 ", hfile, " > ", htempfile))
    system (paste0("gzip --force ", hfile))
    
}

# pbal 
if(file.exists(paste0(pfile,".gz"))) {
    print(paste0("PBALclust result exists ", pfile, ".gz"))
    system (paste0 ("gunzip ", pfile,".gz"))
    ptempfile <- paste0(outdir,"/PBALclust_temp_maxk_",Max_K,".tsv")
    system (paste0 ("cut -f1 ", pfile, " > ", idtempfile))    
    system (paste0 ("cut -f2-11 ", pfile, " > ", ptempfile))
    system (paste0("gzip --force ", pfile))  
}
# For now I am not using Pearsonclust as initialization
# densitycut 
if(file.exists(paste0(dfile,".gz"))) {
    print(paste0("densitycut result exists ", dfile, ".gz"))
    system (paste0 ("gunzip ", dfile,".gz"))
    dtempfile <- paste0(outdir,"/densitycut_temp_maxpc_",maxpc,".tsv")
    system (paste0 ("cut -f1 ", dfile, " > ", idtempfile))    
    system (paste0 ("cut -f2 ", dfile, " > ", dtempfile))
    system (paste0("gzip --force ", dfile))      
} 
 
 # I should add the name of PBAL or HCLUST - added
# Write the file only when at least one of the files exists
if(file.exists(paste0(hfile,".gz")) || file.exists(paste0(pfile,".gz")) || file.exists(paste0(dfile,".gz")))  {
    command <- paste0 ("paste ", idtempfile, " ", htempfile, " ", ptempfile, " ", dtempfile, " > ", outfile)
    system(command)
    system (paste0("gzip --force ", outfile))    
}    

if (ptempfile != "")  system (paste0("rm ", ptempfile))
if (htempfile != "")  system (paste0("rm ", htempfile))
if (dtempfile != "")  system (paste0("rm ", dtempfile))  
if (idtempfile != "")  system (paste0("rm ", idtempfile))

# PYTHON3 <- "/home/mandronescu/.local/centos6/anaconda3/bin/python3"

# If there is a true clusters file, then run evaluation software
if (!is.null(true_clusters_file)) {

    if (hclust_region_crash == 0 && hclust_region_bestpartition_crash == 0) {
        print("Calling evaluation software for Hclust")
        command <- paste0("python3 ", eval_soft, " --true_clusters_file ", true_clusters_file, " --true_prevalences ", args$true_prevalences, " --predicted_clusters_file ", hfile, ".gz --clusters_are_probabilities False --results_file ", outdir, "/results_Hclust.txt")
        print(command)
        system(command)
    }
    
    if (PBAL_crash ==0 && PBALclust_bestpartition_crash == 0) {
        print("Calling evaluation software for PBALclust")
        command <- paste0("python3 ", eval_soft, " --true_clusters_file ", true_clusters_file, " --true_prevalences ", args$true_prevalences, " --predicted_clusters_file ", pfile, ".gz --clusters_are_probabilities False --results_file ", outdir, "/results_PBALclust.txt")
        print(command)
        system(command)
    }    

    if (Pearson_crash ==0 && Pearsonclust_bestpartition_crash == 0) {
        print("Calling evaluation software for Pearsonclust")
        command <- paste0("python3 ", eval_soft, " --true_clusters_file ", true_clusters_file, " --true_prevalences ", args$true_prevalences, " --predicted_clusters_file ", pfile, ".gz --clusters_are_probabilities False --results_file ", outdir, "/results_Pearsonclust.txt")
        print(command)
        system(command)
    }    


    if(file.exists(paste0(dfile,".gz"))) {
        print("Calling evaluation software for densityCut")
        command <- paste0("python3 ", eval_soft, " --true_clusters_file ", true_clusters_file, " --true_prevalences ", args$true_prevalences, " --predicted_clusters_file ", dfile, ".gz --clusters_are_probabilities False --results_file ", outdir, "/results_densitycut.txt")
        print(command)
        system(command)  
    }       
}


