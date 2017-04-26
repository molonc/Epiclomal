
#======================
# libraries
#======================

suppressMessages(library(argparse))

### libraries to find the best number of clusters
# suppressMessages(library(factoextra))
# suppressMessages(library(NbClust))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory") 
parser$add_argument("--methylation_file", type="character", help="Path to methylation data") 
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--max_k", type="integer",default=5, help="maximum number of clusters to be considered when cutting the tree") 

args <- parser$parse_args() 

print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)
input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file


# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
 
# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci5000_clones3_cells100_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.85_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci5000_clones3_cells100_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.85_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci5000_clones3_cells100_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.85_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci5000_clones3_cells100_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.85_gpbrandom_dirpar1_1_nregs50_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells50_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.95_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells50_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.95_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells50_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.95_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells50_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.95_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"


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

#======================
# hiearchical clustering considering Euclidean distances and complete linkage 
#======================

  if (R == 1){

    print("One region, CpG based hiearchical clustering")

    hcluster <- hclust(dist(input_CpG_data,method="euclidean"),method = "complete")
    
    # defining some clusters
    mycl <- cutree(hcluster, k=1:args$max_k)
    
    possible_clusters <- cbind(rownames(input_CpG_data),mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id",paste0("num_clusters_",1:args$max_k))
  
    ofile <- paste0(outdir,"/hclust_clusters_CpG_based_maxk_",args$max_k,".tsv") 
    write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    system(paste0("gzip --force ", ofile))
    
    ofile <- paste0(outdir,"/hclust_cell_order_CpG_based_maxk_",args$max_k,".tsv")
    write.table(hcluster$order, file=ofile , sep="\t", col.names=FALSE, quote=FALSE)
    system(paste0("gzip --force ", ofile))
        
    #save(hcluster, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_hclust_R_object_CpG_based.Rdata"))

   }

  if (R > 1){

    print("More than one region, CpG based hiearchical clustering")
    
    # ### if 95% of the data is missing there is a high chance that for some pairs of cells there won't be any CpG site with data on both cells 
    # ### and, therefore, we won't be able to calculate obtain a hierarchical clustering as the distance/dissimilarity matrix has to be complete
    # 
    # ### verifying how many pairs of cells don't have any CpG site with data on both cells
    # 
    # commom.CpGs.pair <- function(v1,v2){
    #   na.idx <- is.na(v1) | is.na(v2) 
    #   v1a  <- v1[!na.idx]
    #   v2a  <- v2[!na.idx]
    #   l.na <- (sum(!na.idx)) ## = length(v1a) = length(v2a), the number of entries with data on both vectors
    #   return(l.na)
    #   
    # }
    # 
    # commom.CpGs.all <- function(d){ ### d is matrix where the rows correspond to cells and columns to CpGs
    #   dist.data <- matrix(NA,nrow=nrow(d),ncol=nrow(d))
    #   rownames(dist.data) <- rownames(d)
    #   colnames(dist.data) <- rownames(d)
    #   for (i in 1:nrow(d)){
    #     for(j in 1:nrow(d)){
    #       dist.data[i,j] <- commom.CpGs.pair(v1=d[i,],v2=d[j,])
    #     }
    #   }
    #   return(dist.data)
    # }
    # 
    # tmp = commom.CpGs.all(d=input_CpG_data)
    # sum(tmp == 0)
    # rm(tmp)
    # ### In one example when 95% of the data is missing 8% of cell pairs had no CpG with data on both of them
    
    ### there is also the option of doing as Tony did: calculating a dist matrix and then applying hclust on dist of dist,
    ### which is still something I didn't know people do, but it seems they do...
    # tmp = dist(input_CpG_data,method="euclidean")
    # tmp = as.matrix(tmp)
    # plot(hclust(dist(tmp)))
    
    hcluster <- hclust(dist(input_CpG_data,method="euclidean"),method = "complete")
    
    # ### finding optimal number of clusters
    
    # ### example using the iris data set
    # data<-iris[,-c(5)]
    # diss_matrix<- dist(data, method = "euclidean", diag=FALSE)
    # NbClust(data, diss=diss_matrix, distance = NULL, min.nc=2, max.nc=6, 
    #         method = "ward.D2", index = "all") 
    # 
    # ### not working for our data set, not sure why 
    # # diss_matrix <- dist(input_CpG_data,method="euclidean",diag=FALSE)
    # # NbClust(t(input_CpG_data), diss = diss_matrix,distance=NULL, min.nc=2, max.nc=6,,method = "complete",index = "all") 
    # 
    # fviz_nbclust(input_CpG_data, hcut, method = "wss") +
    #   geom_vline(xintercept = 3, linetype = 2)
    
    
    # defining some clusters
    mycl <- cutree(hcluster, k=1:args$max_k)
    
    possible_clusters <- cbind(rownames(input_CpG_data),mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id",paste0("num_clusters_",1:args$max_k))
    
    ofile <- paste0(outdir,"/hclust_clusters_CpG_based_maxk_",args$max_k,".tsv") 
    write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    system(paste0("gzip --force ", ofile))    
    
    ofile <- paste0(outdir,"/hclust_cell_order_CpG_based_maxk_",args$max_k,".tsv") 
    write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste0("gzip --force ", ofile))    
    
    # save(hcluster, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_hclust_R_object_CpG_based.Rdata"))

    print("More than one region, region based hiearchical clustering")
    
    #======================
    # extracting the mean methylation for each region in each cell 
    # this will result in matrix with N cells by R regions - regions are columns and cells the lines 
    #======================
    
    mean_meth_matrix <- t(apply(input_CpG_data,1,extract_mean_meth_per_cell,region_coord=input_regions))

    hcluster <- hclust(dist(mean_meth_matrix ,method="euclidean"),method = "complete")
    
    ### finding the best number of clusters
    # ### working in this case 
    # diss_matrix <- dist(mean_meth_matrix,method="euclidean",diag=FALSE)
    # NbClust(mean_meth_matrix, diss = diss_matrix,distance=NULL, min.nc=2, max.nc=6,,method = "complete",index = "all")
    # 
    # fviz_nbclust(input_CpG_data, hcut, method = "wss") +
    #   geom_vline(xintercept = 3, linetype = 2)
    
    # defining some clusters
    mycl <- cutree(hcluster, k=1:args$max_k)
    
    possible_clusters <- cbind(rownames(input_CpG_data),mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id",paste0("num_clusters_",1:args$max_k))
    
    ofile <- paste0(outdir,"/hclust_clusters_region_based_maxk_",args$max_k,".tsv") 
    write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    system(paste0("gzip --force ", ofile))
    
    ofile <- paste0(outdir,"/hclust_cell_order_region_based_maxk_",args$max_k,".tsv")
    write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste0("gzip --force ", ofile))    
    
    # save(hcluster, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_hclust_R_object_region_based.Rdata"))
  
  }   

### Tony's (PBAL manuscript) clustering approach 

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

hcluster_T <- hclust(dist(dist.PBAL(d=input_CpG_data),method="euclidean"),method = "ward.D2")

# ### Choosing number of clusters
# ### working in this case
# diss_matrix <- dist(dist.PBAL(d=input_CpG_data),method="euclidean",diag=FALSE)
# NbClust(dist.PBAL(d=input_CpG_data), diss = diss_matrix,distance=NULL, min.nc=2, max.nc=6,,method = "complete",index = "ch")
# fviz_nbclust(input_CpG_data, hcut, method = "wss") +
# geom_vline(xintercept = 3, linetype = 2)

# defining some clusters
mycl_T <- cutree(hcluster_T, k=1:args$max_k)

possible_clusters_T <- cbind(rownames(input_CpG_data),mycl_T)
possible_clusters_T <- as.data.frame(possible_clusters_T)
colnames(possible_clusters_T) <- c("cell_id",paste0("num_clusters_",1:args$max_k))

ofile <- paste0(outdir,"/PBALclust_clusters_CpG_based_maxk_",args$max_k,".tsv") 
write.table(possible_clusters_T, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
system(paste0("gzip --force ", ofile))
    
ofile <- paste0(outdir,"/PBALclust_cell_order_CpG_based_maxk_",args$max_k,".tsv")
write.table(hcluster_T$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
system(paste0("gzip --force ", ofile))

#save(hcluster_T, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_PBALclust_R_object_CpG_based.Rdata"))




