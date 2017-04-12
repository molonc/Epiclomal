
#======================
# libraries
#======================

suppressMessages(library(argparse))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--path_directory", type="character", help="Path to working directory") 
parser$add_argument("--methylation_file", type="character", help="Path to methylation data") 
parser$add_argument("--max_k", type="integer",default=5, help="maximum number of clusters to be considered when cutting the tree") 
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")

args <- parser$parse_args() 

print(args)

input_CpG_data_file <- paste0(args$path_directory,args$methylation_file)
input_regions_file <- paste0(args$path_directory,args$regions_file)


# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
 
# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.7_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"


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
  
    write.table(possible_clusters, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_maxk_",args$max_k,"_hclust_clusters_CpG_based.tsv") , sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    write.table(hcluster$order, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_hclust_cell_order_CpG_based.tsv") , sep="\t", col.names=FALSE, quote=FALSE)
    save(hcluster, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_hclust_R_object_CpG_based.Rdata"))

   }

  if (R > 1){

    print("More than one region, CpG based hiearchical clustering")
    
    hcluster <- hclust(dist(input_CpG_data,method="euclidean"),method = "complete")
    
    # defining some clusters
    mycl <- cutree(hcluster, k=1:args$max_k)
    
    possible_clusters <- cbind(rownames(input_CpG_data),mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id",paste0("num_clusters_",1:args$max_k))
    
    write.table(possible_clusters, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_maxk_",args$max_k,"_hclust_clusters_CpG_based.tsv") , sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    write.table(hcluster$order, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_hclust_cell_order_CpG_based.tsv") , sep="\t", col.names=FALSE, quote=FALSE)
    save(hcluster, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_hclust_R_object_CpG_based.Rdata"))

    print("More than one region, region based hiearchical clustering")
    
    #======================
    # extracting the mean methylation for each region in each cell 
    # this will result in matrix with N cells by R regions - regions are columns and cells the lines 
    #======================
    
    mean_meth_matrix <- t(apply(input_CpG_data,1,extract_mean_meth_per_cell,region_coord=input_regions))

    hcluster <- hclust(dist(mean_meth_matrix ,method="euclidean"),method = "complete")
    
    # defining some clusters
    mycl <- cutree(hcluster, k=1:args$max_k)
    
    possible_clusters <- cbind(rownames(input_CpG_data),mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id",paste0("num_clusters_",1:args$max_k))
    
    write.table(possible_clusters, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_maxk_",args$max_k,"_hclust_clusters_region_based.tsv") , sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    write.table(hcluster$order, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_hclust_cell_order_region_based.tsv") , sep="\t", col.names=FALSE, quote=FALSE)
    save(hcluster, file=paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_hclust_R_object_region_based.Rdata"))
  
         }   



