
#======================
# libraries
#======================

### for testing I used tehse paths for the libraries and R in /gsc/software/linux-x86_64-centos6/R-3.3.2/bin/R
# library(pcaMethods,lib.loc="/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.2")
# library(densitycut,lib.loc = "/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.3/" )

suppressMessages(library(argparse))
suppressMessages(library(Rcpp))
suppressMessages(library(pcaMethods))
suppressMessages(library(densitycut))
suppressMessages(library(pheatmap))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--methylation_file", type="character", help="Path to methylation data")
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--max_PC", type="integer",default=20, help="maximum number of PC components")
parser$add_argument("--maxit", type="integer",default=100, help="maximum number iterations")
# MA: 30Jan 2019: adding an optional imputation step
parser$add_argument("--impute", default="0", type="integer",help="If it is 1, impute with the average per region/locus, if it is 0 do nothing.")

args <- parser$parse_args()

print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)
input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file
impute <- args$impute

# input_CpG_data_file <- "data_incomplete.tsv.gz"
# input_regions_file <- "regions_file.tsv.gz"

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# #outdir <- "/Users/cdesouza/Documents/synthetic_data_old/output_loci1000_clones3_cells100_prev0.2_0.5_0.3_errpb0.001_0.001_mispb0_gpbrandom_dirpar1_1_nregs4_rsize-equal_rnonequal-uniform/data"
#
#  input_CpG_data_file <- "/genesis/shahlab/csouza/BS-seq/whole_genome_single_cell/synthetic_data_tests/data_old_way/data/data_incomplete.tsv.gz"
#  input_regions_file <- "/genesis/shahlab/csouza/BS-seq/whole_genome_single_cell/synthetic_data_tests/data_old_way/data/regions_file.tsv.gz"
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
cached_data <- gsub(".tsv.gz", ".RData.gz", input_CpG_data_file)
if (file.exists(cached_data)) {
  print("loading cached data")
  load(cached_data)
} else {
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

  mean_meth_matrix <- t(apply(input_CpG_data,1,extract_mean_meth_per_cell,region_coord=input_regions))

  save(input_CpG_data, input_regions, mean_meth_matrix, file = cached_data, compress = "gzip")
}


R <- dim(input_regions)[1] ## number of regions
M <- dim(input_CpG_data)[2] ## number of loci

print("number of regions:")
print(R)
print("number of loci:")
print(M)

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

  print("More than one region, region based densitycut")

  #======================
  # extracting the mean methylation for each region in each cell
  # this will result in matrix with N cells by R regions - regions are columns and cells the lines
  #======================

  # this code is redundant with the code in hclust.R, TO FIX
  if (impute == 1) {
    imputed_file <- paste0(outdir,"/region_based_imputed.RData.gz")
    if (file.exists(imputed_file)) {
      print ("Reading the imputed file")
      load(imputed_file)
      print (" ... done.")
    } else {
      # to remove
      #input <- input[1:100,1:5]

      # replace with average values, for each col
      # get directory of current script
      initial.options <- commandArgs(trailingOnly = FALSE)
      file.arg.name <- "--file="
      script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
      script.basename <- dirname(script.name)
      sourceCpp(paste(sep="/", script.basename, "impute.cpp"))
      print("Per region, replacing NAs with average values")
      mean_meth_matrix <- impute_means(mean_meth_matrix)
      print(" ... done.")

      # eliminate the empty rows (features)
      mean_meth_matrix <- mean_meth_matrix[ rowSums(mean_meth_matrix)!=0, ]
      save(mean_meth_matrix, file = imputed_file, compress = "gzip")
    }
  }




  max_comp <- min(args$max_PC,R)

  print("number of PC components:")
  print(max_comp)

  t <- try(pca( mean_meth_matrix ,method="nipals",nPcs=max_comp))
  if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
    print("Stop! At least one cell has NO data across all regions, can't do PCA")
    stop()
    } else {
    pc <- t }

  pc_scores <- scores(pc)

  #print(head(pc_scores))
  #plot(pc_scores[,1],pc_scores[,2],col=true_membership)

  cluster.out <- DensityCut(pc_scores,maxit=args$maxit) # DensityCut clustering analysis CS: increased max number of iterations from 50 to 100
  checking_warning <- capture.output(cluster.out)

  #print(checking_warning)

  if(checking_warning[1] == "WARNING! not converged "){
    print("densitycut didn't converge, not saving results!")
    } else {

    print("densitycut converged!")


    #print(cluster.out$cluster)

    #col <- AssignLabelColor(label=cluster.out$cluster, col=distinct.col) # Assign colour to clusters
    #NeatPlot(x=pc_scores, col=col) # Scatter plots

    # possible_clusters <- cbind(as.numeric(rownames(input_CpG_data)),cluster.out$cluster)
    # MA: removing as.numeric, otherwise it replaces the cell ids with NAs whenever the cell ids are not numeric

    possible_clusters <- cbind(rownames(input_CpG_data),cluster.out$cluster)
    colnames(possible_clusters) <- c("cell_id","DensityCut")

    #print(possible_clusters)

    ### Plotting DensityCut clusters
    print("Plotting DensityCut clusters")

    inf_clusters_order <- order(as.integer(cluster.out$cluster))

    #print(inf_clusters_order)

    annotation_row <- as.data.frame(cluster.out$cluster)
    colnames(annotation_row) <- "inferred clusters"
    annotation_row$`inferred clusters` <- as.factor(annotation_row$`inferred clusters`)

    rownames(annotation_row) <- rownames(mean_meth_matrix) ### CS: I added this and the bug regarding --> Error in check.length("fill") :
                                                         ### 'gpar' element 'fill' must not be length 0
                                                         ### Calls: pheatmap ... rectGrob -> grob -> gpar -> validGP -> check.length
                                                         ### is now fixed

    #print(annotation_row)
    #print(str(annotation_row))

    if(!file.exists(paste0(outdir,"/DensityCut_PLOT.pdf"))){
      pheatmap(mean_meth_matrix[inf_clusters_order,],cluster_rows = FALSE,cluster_cols=FALSE,
                cellwidth = 8,
                annotation_row = annotation_row,
                cellheight = 8,fontsize = 8,
                #clustering_distance_rows = "euclidean",
                #clustering_method = "ward.D2",
                main = paste0("DensityCut"),
                show_colnames=FALSE,
                annotation_names_row = FALSE,
                filename = paste0(outdir,"/DensityCut_PLOT.pdf"))
    }

  ### end of plotting

  ofile <- paste0(outdir,"/DensityCut_clusters_Region_based_maxPC_",max_comp,".tsv")
  write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
  system(paste0("gzip --force ", ofile))


  }

}

print("DONE!")




