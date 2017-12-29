
#======================
# libraries
#======================

### for testing I used this paths for the libraries and R in /gsc/software/linux-x86_64-centos6/R-3.3.2/bin/R
# library(pcaMethods,lib.loc="/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.2")
# library(densitycut,lib.loc = "/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.3/" )

suppressMessages(library(argparse))
suppressMessages(library(pcaMethods))
suppressMessages(library(densitycut))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory") 
parser$add_argument("--methylation_file", type="character", help="Path to methylation data") 
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--max_PC", type="integer",default=20, help="maximum number of PC components") 

args <- parser$parse_args() 

print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)
input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file

# input_CpG_data_file <- "data_incomplete.tsv.gz"
# input_regions_file <- "regions_file.tsv.gz"

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# # outdir <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data"
# 
#  input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform_seed_2/data/data_incomplete.tsv"
#  input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform_seed_2/data/regions_file.tsv"
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

print(R)
print(M)

#======================
# hiearchical clustering considering Euclidean distances and complete linkage 
#======================

if (R == 1){
  
  print("One region, cannot do Region-based densityCut")
  
  # print("One region, CpG based densityCut")
  # 
  # pc <- pca(input_CpG_data,method="nipals",nPcs=args$max_PC)
  # 
  # pc_scores <- scores(pc)
  # 
  # print(head(pc_scores))
  # 
  # #plot(pc_scores[,1],pc_scores[,2],col=true_membership)
  # 
  # cluster.out <-  DensityCut(pc_scores) # DensityCut clustering analysis
  # 
  # #col <- AssignLabelColor(label=cluster.out$cluster, col=distinct.col) # Assign colour to clusters
  # #NeatPlot(x=pc_scores, col=col) # Scatter plots
  # 
  # print(cluster.out$cluster)
  # 
  # possible_clusters <- cbind(as.numeric(rownames(input_CpG_data)),cluster.out$cluster)
  # colnames(possible_clusters) <- c("cell_id","densitycut")
  # 
  # 
  # ofile <- paste0(outdir,"/density_clusters_CpG_based_maxPC_",args$max_PC,".tsv") 
  # write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
  # system(paste0("gzip --force ", ofile))   
  
}


if (R > 1){
  
  # print("More than one region, CpG based densityCut")
  # 
  # pc <- pca(input_CpG_data,method="nipals",nPcs=args$max_PC)
  # 
  # sum(apply(input_CpG_data,2,function(x){sum(is.na(x))}) == 50)
  # 
  # pc_scores <- scores(pc)
  # 
  # print(head(pc_scores))
  # 
  # #plot(pc_scores[,1],pc_scores[,2],col=true_membership)
  # 
  # cluster.out <-  DensityCut(pc_scores) # DensityCut clustering analysis
  # 
  # #col <- AssignLabelColor(label=cluster.out$cluster, col=distinct.col) # Assign colour to clusters
  # #NeatPlot(x=pc_scores, col=col) # Scatter plots
  # 
  # print(cluster.out$cluster)
  # 
  # possible_clusters <- cbind(as.numeric(rownames(input_CpG_data)),cluster.out$cluster)
  # colnames(possible_clusters) <- c("cell_id","densitycut")
  #   
  # 
  # ofile <- paste0(outdir,"/density_clusters_CpG_based_maxPC_",args$max_PC,".tsv") 
  # write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
  # system(paste0("gzip --force ", ofile))    
  # 
  # rm(cluster.out)
  
  print("More than one region, region based hiearchical clustering")
  
  #======================
  # extracting the mean methylation for each region in each cell 
  # this will result in matrix with N cells by R regions - regions are columns and cells the lines 
  #======================
  
  mean_meth_matrix <- t(apply(input_CpG_data,1,extract_mean_meth_per_cell,region_coord=input_regions))
  
  max_comp <- min(args$max_PC,R)
  
  t <- try(pca( mean_meth_matrix ,method="nipals",nPcs=max_comp))
  if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
    print("Stop! At least one cell has NO data across all regions, can't do PCA")
    stop()
    }else{
    pc <- pca( mean_meth_matrix ,method="nipals",nPcs=max_comp) }
  
  pc_scores <- scores(pc)
  
  #print(head(pc_scores))
  
  #plot(pc_scores[,1],pc_scores[,2],col=true_membership)
  
  checking_warning <- capture.output(DensityCut(pc_scores))
  
  if(checking_warning[1] == "WARNING! not converged "){
    print("densitycut didn't converge, not saving results")}else{
  
  cluster.out <-  DensityCut(pc_scores) # DensityCut clustering analysis
  
  #print(cluster.out$cluster)
  
  #col <- AssignLabelColor(label=cluster.out$cluster, col=distinct.col) # Assign colour to clusters
  #NeatPlot(x=pc_scores, col=col) # Scatter plots
  
  possible_clusters <- cbind(as.numeric(rownames(input_CpG_data)),cluster.out$cluster)
  colnames(possible_clusters) <- c("cell_id","densitycut")
  
  ofile <- paste0(outdir,"/densityCut_clusters_Region_based_maxPC_",max_comp,".tsv") 
  write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
  system(paste0("gzip --force ", ofile))   
  
  }
  
}


