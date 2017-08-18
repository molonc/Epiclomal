
#======================
# libraries
#======================

suppressMessages(library(argparse))

### libraries to find the best number of clusters
#suppressMessages(library(factoextra))
suppressMessages(library(NbClust))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory") 
parser$add_argument("--methylation_file", type="character", help="Path to methylation data") 
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--max_k", type="integer",default=5, help="maximum number of clusters to be considered when cutting the tree") 
#parser$add_argument("--include_true_clone_membership",default="no", type="character",help="yes or no; for analysis of real data it should be no")
parser$add_argument("--true_clone_membership_file", type="character",help="Path to true clone membership file")

args <- parser$parse_args() 

print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)
input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file
true_clusters_file <- args$true_clone_membership_file

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# outdir <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data"

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform_seed_2/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform_seed_2/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform_seed_2/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform_seed_2/data/true_clone_membership.tsv"

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
    
    if  (!is.null(args$true_clone_membership_file)){
      
      # true clone cell membership
      tmp <- read.csv(true_clusters_file,sep="\t",header=TRUE,check.names=FALSE)
      true_membership <- as.matrix(tmp[,-1])
      rm(tmp)  
      
      possible_clusters <- cbind(rownames(input_CpG_data),true_membership,mycl)
      possible_clusters <- as.data.frame(possible_clusters)
      colnames(possible_clusters) <- c("cell_id","true_membership",paste0("num_clusters_",1:Max_K))
    } else{
      possible_clusters <- cbind(rownames(input_CpG_data),mycl)
      possible_clusters <- as.data.frame(possible_clusters)
      colnames(possible_clusters) <- c("cell_id",paste0("num_clusters_",1:Max_K))
    }
    
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
  
  # print("More than one region, CpG based hiearchical clustering")
  # 
  # # ### if 95% of the data is missing there is a high chance that for some pairs of cells there won't be any CpG site with data on both cells 
  # # ### and, therefore, we won't be able to calculate obtain a hierarchical clustering as the distance/dissimilarity matrix has to be complete
  # # 
  # # ### verifying how many pairs of cells don't have any CpG site with data on both cells
  # # 
  # # commom.CpGs.pair <- function(v1,v2){
  # #   na.idx <- is.na(v1) | is.na(v2) 
  # #   v1a  <- v1[!na.idx]
  # #   v2a  <- v2[!na.idx]
  # #   l.na <- (sum(!na.idx)) ## = length(v1a) = length(v2a), the number of entries with data on both vectors
  # #   return(l.na)
  # #   
  # # }
  # # 
  # # commom.CpGs.all <- function(d){ ### d is matrix where the rows correspond to cells and columns to CpGs
  # #   dist.data <- matrix(NA,nrow=nrow(d),ncol=nrow(d))
  # #   rownames(dist.data) <- rownames(d)
  # #   colnames(dist.data) <- rownames(d)
  # #   for (i in 1:nrow(d)){
  # #     for(j in 1:nrow(d)){
  # #       dist.data[i,j] <- commom.CpGs.pair(v1=d[i,],v2=d[j,])
  # #     }
  # #   }
  # #   return(dist.data)
  # # }
  # # 
  # # tmp = commom.CpGs.all(d=input_CpG_data)
  # # sum(tmp == 0)
  # # rm(tmp)
  # # ### In one example when 95% of the data is missing 8% of cell pairs had no CpG with data on both of them
  # 
  # ### there is also the option of doing as Tony did: calculating a dist matrix and then applying hclust on dist of dist,
  # ### which is still something I didn't know people do, but it seems they do...
  # # tmp = dist(input_CpG_data,method="euclidean")
  # # tmp = as.matrix(tmp)
  # # plot(hclust(dist(tmp)))
  # 
  # pairwisedist <- dist(input_CpG_data,method="euclidean")
  # print(sum(is.na(pairwisedist) == TRUE))
  # 
  # if(sum(is.na(pairwisedist)==TRUE) == 0){
  #   
  #   hclust_CpG_crash <- 0
  #   write.table(hclust_CpG_crash,file=paste0(outdir,"/hclust_CpG_crash.tsv"),row.names=FALSE,col.names=FALSE)
  #   
  #   hcluster <- hclust(pairwisedist,method = "complete")
  #   
  #   mycl <- cutree(hcluster, k=1:Max_K)
  #   
  #   if  (!is.null(args$true_clone_membership_file)){
  #     
  #     # true clone cell membership
  #     tmp <- read.csv(true_clusters_file,sep="\t",header=TRUE,check.names=FALSE)
  #     true_membership <- as.matrix(tmp[,-1])
  #     rm(tmp)  
  #     
  #     possible_clusters <- cbind(rownames(input_CpG_data),true_membership,mycl)
  #     possible_clusters <- as.data.frame(possible_clusters)
  #     colnames(possible_clusters) <- c("cell_id","true_membership",paste0("num_clusters_",1:Max_K))
  #   } else{
  #     possible_clusters <- cbind(rownames(input_CpG_data),mycl)
  #     possible_clusters <- as.data.frame(possible_clusters)
  #     colnames(possible_clusters) <- c("cell_id",paste0("num_clusters_",1:Max_K))
  #   }
  #   
  #   
  #   ### using package clusterCrit
  #   #cl <- kmeans(input_CpG_data,3)
  #   #intCriteria(traj=input_CpG_data,part=as.integer(cl$cluster),crit="Silhouette") ### it doesn't work because the entries of input_CpG_data are zeros and ones
  #   
  #   t <- try(NbClust(input_CpG_data, diss = pairwisedist,distance=NULL, min.nc=2, max.nc=Max_K,method = "complete",index = "cindex"))
  #   if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
  #     print("can't find a best partition")
  #     error_ch_index <- 1 }
  #   else {
  #     error_ch_index <- 0
  #     hcluster_Nb <- NbClust(input_CpG_data, diss = pairwisedist,distance=NULL, min.nc=2, max.nc=Max_K,method = "complete",index = "cindex")
  #     print(hcluster_Nb)}
  #   if(error_ch_index == 1){
  #     write.table(error_ch_index,file=paste0(outdir,"/hclust_CpGbased_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)
  #     
  #     ofile <- paste0(outdir,"/hclust_clusters_CpG_based_maxk_",Max_K,".tsv")
  #     write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
  #     system(paste0("gzip --force ", ofile))
  #     
  #     ofile <- paste0(outdir,"/hclust_cell_order_CpG_based_maxk_",Max_K,".tsv")
  #     write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
  #     system(paste0("gzip --force ", ofile))
  #   }
  #   
  #   if(error_ch_index == 0){
  #     
  #     write.table(error_ch_index,file=paste0(outdir,"/hclust_CpGbased_bestpartition_crash.tsv"),row.names=FALSE,col.names=FALSE)
  #     
  #     best_cluster <- hcluster_Nb$Best.partition
  #     
  #     possible_clusters <- cbind(possible_clusters,best_cluster)
  #     colnames(possible_clusters) <- c(colnames(possible_clusters)[1:(dim(possible_clusters)[2]-1)],paste0("best_cluster_",hcluster_Nb$Best.nc[1]))
  #     
  #     ofile <- paste0(outdir,"/hclust_clusters_CpG_based_maxk_",Max_K,".tsv") 
  #     write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
  #     system(paste0("gzip --force ", ofile))    
  #     
  #     ofile <- paste0(outdir,"/hclust_cell_order_CpG_based_maxk_",Max_K,".tsv") 
  #     write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
  #     system(paste0("gzip --force ", ofile))    
  #     
  #     rm(hcluster_Nb)
  #     
  #   }
  #   
  #   rm(hcluster)
  #   
  # }
  # 
  # if(sum(is.na(pairwisedist)==TRUE) > 0){
  #   print("some pairs of cells have no CpG with data in common")
  #   
  #   hclust_CpG_crash <- 1
  #   write.table(hclust_CpG_crash,file=paste0(outdir,"/hclust_CpG_crash.tsv"),row.names=FALSE,col.names=FALSE)
  # }
  # 
  
  print("More than one region, region based hiearchical clustering")
  
  #======================
  # extracting the mean methylation for each region in each cell 
  # this will result in matrix with N cells by R regions - regions are columns and cells the lines 
  #======================
  
  mean_meth_matrix <- t(apply(input_CpG_data,1,extract_mean_meth_per_cell,region_coord=input_regions))
  
  pairwisedist_region <- dist(mean_meth_matrix ,method="euclidean")
  
  if(sum(is.na(pairwisedist_region)==TRUE) == 0){
    
    hclust_region_crash <- 0
    write.table(hclust_region_crash,file=paste0(outdir,"/hclust_region_crash.tsv"),row.names=FALSE,col.names=FALSE)
    
    hcluster <- hclust(pairwisedist_region,method = "complete")
    
    # defining some clusters
    mycl <- cutree(hcluster, k=1:Max_K)
    
    if  (!is.null(args$true_clone_membership_file)){
      
      # true clone cell membership
      tmp <- read.csv(true_clusters_file,sep="\t",header=TRUE,check.names=FALSE)
      true_membership <- as.matrix(tmp[,-1])
      rm(tmp)  
      
      possible_clusters <- cbind(rownames(input_CpG_data),true_membership,mycl)
      possible_clusters <- as.data.frame(possible_clusters)
      colnames(possible_clusters) <- c("cell_id","true_membership",paste0("num_clusters_",1:Max_K))
    } else{
      possible_clusters <- cbind(rownames(input_CpG_data),mycl)
      possible_clusters <- as.data.frame(possible_clusters)
      colnames(possible_clusters) <- c("cell_id",paste0("num_clusters_",1:Max_K))
    }
    
    ### finding the best number of clusters
    # ### working in this case 
    #diss_matrix <- dist(mean_meth_matrix,method="euclidean",diag=FALSE)
    
    ## t <- try(NbClust(mean_meth_matrix, diss = diss_matrix,distance=NULL, min.nc=2, max.nc=Max_K,method = "complete",index = "ch"))
    t <- try(NbClust(mean_meth_matrix, diss = pairwisedist_region, distance=NULL, min.nc=2, max.nc=Max_K,method = "complete",index = "ch")) ### changing to cindex as cindex also works for CpG based clustering
    if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
      print("can't use ch index") 
      error_ch_index <- 1 }
    else { 
      error_ch_index <- 0
      hcluster_Nb <- NbClust(mean_meth_matrix, diss = pairwisedist_region,distance=NULL, min.nc=2, max.nc=Max_K,method = "complete",index = "ch")
    print(hcluster_Nb)}
    
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
  
  # defining some clusters
  mycl_T <- cutree(hcluster_T, k=1:Max_K)
  
  if  (!is.null(args$true_clone_membership_file)){
    
    # true clone cell membership
    tmp <- read.csv(true_clusters_file,sep="\t",header=TRUE,check.names=FALSE)
    true_membership <- as.matrix(tmp[,-1])
    rm(tmp)  
    
    possible_clusters_T <- cbind(rownames(input_CpG_data),true_membership,mycl_T)
    possible_clusters_T <- as.data.frame(possible_clusters_T)
    colnames(possible_clusters_T) <- c("cell_id","true_membership",paste0("num_clusters_",1:Max_K))
  } else{
    possible_clusters_T <- cbind(rownames(input_CpG_data),mycl_T)
    possible_clusters_T <- as.data.frame(possible_clusters_T)
    colnames(possible_clusters_T) <- c("cell_id",paste0("num_clusters_",1:Max_K))
  }
  
  ## t <- try(NbClust(mean_meth_matrix, diss = diss_matrix,distance=NULL, min.nc=2, max.nc=Max_K,method = "complete",index = "cindex"))
  t <- try(NbClust(dist_PBAL, diss = diss_matrix_T,distance=NULL, min.nc=2, max.nc=Max_K,method = "ward.D2",index = "ch")) ### changing to cindex as cindex also works for CpG based clustering
  if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
    print("can't find best partition") 
    error_ch_index <- 1 }
  else { 
    error_ch_index <- 0
    hcluster_Nb_T <- NbClust(dist_PBAL, diss = diss_matrix_T,distance=NULL, min.nc=2, max.nc=Max_K,method = "ward.D2",index = "ch")
    print(hcluster_Nb_T)
  }
  
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


# MA: added another file at the end with all the columns from hclust regions and pbal (except the first 2 columns of pbal cell_id and clusters_1
# Note: I am unzipping so I zip again after, I should not zip earlier, TODO
# TODO: skip the best column, this file will be used only for initialization, not for evaluation

hfile <- paste0(outdir,"/hclust_clusters_region_based_maxk_",Max_K,".tsv")
pfile <- paste0(outdir,"/PBALclust_clusters_CpG_based_maxk_",Max_K,".tsv")
outfile <- paste0(outdir, "/hclust_region_PBAL_CpG_clusters_maxk_",Max_K,".tsv")
print (paste0("Hfile ", hfile))
print (paste0("Pfile ", pfile))
print (paste0("Outfile ", outfile))
system (paste0 ("gunzip ", hfile,".gz"))
system (paste0 ("gunzip ", pfile,".gz"))
htempfile <- paste0(outdir,"/hclust_temp_maxk_",Max_K,".tsv")
ptempfile <- paste0(outdir,"/PBALclust_temp_maxk_",Max_K,".tsv")
system (paste0 ("cut -f1-11 ", hfile, " > ", htempfile))
system (paste0 ("cut -f2-11 ", pfile, " > ", ptempfile))
command <- paste0 ("paste ", htempfile, " ", ptempfile, " > ", outfile)
system(command)
system (paste0("rm ", ptempfile))
system (paste0("rm ", htempfile))
system (paste0("gzip --force ", hfile))
system (paste0("gzip --force ", pfile))
system (paste0("gzip --force ", outfile))
