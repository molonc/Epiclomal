# generate output plots for Epiclomal
#   Made each file have its own path, replaced path_dir with out_dir
#   Put all the results in out_dir
#   Added a name for the plot to identify where the clustering is coming from
#   Added order_by_true argument

#======================
# libraries
#======================
library(pheatmap)
suppressMessages(library(argparse))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--copy_number", type="integer", default=0, help="Set to 1 if you want to plot copy number data")
parser$add_argument("--true_clusters", type="integer", default=1, help="Set to 1 if true cluster memberships are available")
parser$add_argument("--out_directory", type="character", help="Path to output directory") 
parser$add_argument("--methylation_file", type="character", help="Path to methylation data") 
parser$add_argument("--copy_number_file", type="character", help="Path to copy number data if available")
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--inferred_clusters_file", type="character", help="Path to inferred clusters. The third column has posterior probabilities.")
# MA: the inferred_clusters_file may have a thurd column with the posterior probabilities
# parser$add_argument("--cell_posterior_probabilities_file", type="character", help="Path to cell posterior probabilities of belonging to a cluster")

parser$add_argument("--true_clusters_file", type="character", help="Path to true clusters if available")
parser$add_argument("--order_by_true", type="integer", default=1, help="If 1, rows are ordered by the true clustering if given, else by the predicted")
parser$add_argument("--name", type="character", help="A name for the final plot")

args <- parser$parse_args() 

print(args)

out_dir <- args$out_directory
dir.create(file.path(out_dir), showWarnings = FALSE)

input_CpG_data_file <- paste0(args$path_directory,args$methylation_file)
input_regions_file <- paste0(args$path_directory,args$regions_file)
inferred_clusters_file <-  paste0(args$path_directory,args$inferred_clusters_file)

#cell_posterior_probabilities_file <-  paste0(args$path_directory,args$cell_posterior_probabilities_file)

if ( args$true_clusters == 1){
  true_clusters_file <- paste0(args$path_directory,args$true_clusters_file) }

if ( args$copy_number == 1){
  input_CN_data_file <- paste0(args$path_directory,args$copy_number_file) }

#  input_CN_data_file <- "/Users/cdesouza/Documents/EPI-91/CN_data_most_variable_CGIs_xeno7_Epiclomal.tsv"
#  input_CpG_data_file <- "/Users/cdesouza/Documents/EPI-91/most_var_CGIs_all_cells_input_Epiclomal_hg19_xeno7.tsv"
#  input_regions_file <- "/Users/cdesouza/Documents/EPI-91/most_var_CGIs_regionIDs_input_Epiclomal_hg19_xeno7.tsv"
#  inferred_clusters_file <- "/Users/cdesouza/Documents/EPI-91/result_Basic_CN_Epiclomal_100repeats_MAXK7_64.tsv"
#  cell_posterior_probabilities_file <- "/Users/cdesouza/Documents/EPI-91/result_Basic_CN_Epiclomal_100repeats_MAXK7_64_posterior_fake.tsv"

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs5_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"

# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs1_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs1_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs1_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci100_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs1_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# 
# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci500_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs1_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci500_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs1_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci500_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs1_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci500_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs1_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# 
# input_CpG_data_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci5000_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs100_rsize-equal_rnonequal-uniform/data/data_incomplete.tsv"
# input_regions_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci5000_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs100_rsize-equal_rnonequal-uniform/data/regions_file.tsv"
# inferred_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci5000_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs100_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"
# true_clusters_file <- "/Users/cdesouza/Documents/synthetic_data/output_loci5000_clones3_cells40_prev0.2_0.5_0.3_errpb0.01_0.01_mispb0.25_gpbrandom_dirpar1_1_nregs100_rsize-equal_rnonequal-uniform/data/true_clone_membership.tsv"


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

# Copy number data
if (args$copy_number == 1 ){
  tmp <- read.csv(input_CN_data_file,sep="\t",header=TRUE,check.names=FALSE)
  input_CN_data <- as.matrix(tmp[,-1])
  rownames(input_CN_data) <- tmp$cell_id
  colnames(input_CN_data) <- sapply(strsplit(colnames(input_CN_data),":"),function(x){return(x[1])}) 
  rm(tmp)
}

# Region coordinates
tmp <- read.csv(input_regions_file,sep="\t",header=TRUE,check.names=FALSE)
input_regions <- as.matrix(tmp[,-1]) + 1 ## adding 1 to match R indexing - previously coordinates were for python starting on zero
colnames(input_regions) <- c("start","end") ## input_regions gives already the columns in input_CpG_data that correspond to which regions considered in the construction of input_CpG_data
rownames(input_regions) <- tmp$region_id
rm(tmp)

# inferred clusters from Epliclomal
tmp <- read.csv(inferred_clusters_file,sep="\t",header=TRUE,check.names=FALSE)
colnames(tmp) <- c("cell_id","inferred_clusters") 
inferred_cell_clusters <- as.matrix(tmp[,2])
rownames(inferred_cell_clusters) <- tmp$cell_id
inferred_cell_clusters <- as.data.frame(inferred_cell_clusters)
colnames(inferred_cell_clusters) <- "inferred_clusters"
#tmp[47:49,] <- tmp[c(48,49,47),]
if(sum(rownames(input_CpG_data) != rownames(inferred_cell_clusters)) > 0){
  stop("order of cell IDs doesn't match")
}
rm(tmp)


# # cell posterior probabilities
# tmp <- read.csv(cell_posterior_probabilities_file,sep="\t",header=TRUE,check.names=FALSE)
# cell_posterior_probabilities <- as.matrix(apply(tmp[,-1],1,max)) ### check the format of this file
# colnames(tmp) <- c("cell_id","inferred_clusters") 
# rownames(cell_posterior_probabilities) <- tmp$cell_id
# cell_posterior_probabilities <- as.data.frame(cell_posterior_probabilities)
# colnames(cell_posterior_probabilities) <- "cell_posteriors"
# if(sum(rownames(input_CpG_data) != rownames(cell_posterior_probabilities)) > 0){
#   stop("order of cell IDs doesn't match")
# }
# rm(tmp)


if ( args$true_clusters == 1){
  tmp <- read.csv(true_clusters_file,sep="\t",header=TRUE,check.names=FALSE)
  colnames(tmp) <- c("cell_id","true_clusters") 
  true_cell_clusters <- as.matrix(tmp[,2])
  rownames(true_cell_clusters) <- tmp$cell_id
  true_cell_clusters <- as.data.frame(true_cell_clusters)
  colnames(true_cell_clusters) <- "true_clusters"
  #tmp[47:49,] <- tmp[c(48,49,47),]
  if(sum(rownames(input_CpG_data) != rownames(true_cell_clusters)) > 0){
    stop("order of cell IDs doesn't match")
  }
  rm(tmp)
}

#======================
# extracting the mean methylation for each region in each cell 
# this will result in matrix with N cells by R regions - regions are columns and cells the lines 
#======================

mean_meth_matrix <- t(apply(input_CpG_data,1,extract_mean_meth_per_cell,region_coord=input_regions))
## checking row and column names
#rownames(mean_meth_matrix)
#colnames(mean_meth_matrix)

#======================
# pheatmap plots
#======================

if(args$true_clusters == 1 && args$order_by_true == 1){
  index <- 1:dim(true_cell_clusters)[1]
  index_gaps <- index[!duplicated(true_cell_clusters[order(true_cell_clusters),])] - 1 
  index_gaps <- index_gaps[which(index_gaps != 0)]  
} else{
  index <- 1:dim(inferred_cell_clusters)[1]
  index_gaps <- index[!duplicated(inferred_cell_clusters[order(inferred_cell_clusters),])] - 1 
  index_gaps <- index_gaps[which(index_gaps != 0)]
}

## annotating the rows by clusters

# MA: we want to add the true clusters annotation no matter whether we order by true or predicted
if(args$true_clusters == 1){
  annotation_row <- cbind(inferred_cell_clusters,true_cell_clusters$true_clusters)
  colnames(annotation_row) <- c("inferred clusters", "true clusters")  
  annotation_row$`inferred clusters` <- as.factor(annotation_row$`inferred clusters`)
  annotation_row$`true clusters` <- as.factor(annotation_row$`true clusters`)
}else{
  annotation_row <- inferred_cell_clusters
  colnames(annotation_row) <-"inferred clusters"
  annotation_row$`inferred clusters` <- as.factor(annotation_row$`inferred clusters`)
  
  # ### adding posterior probabilities
  # annotation_row <- cbind(inferred_cell_clusters,cell_posterior_probabilities$cell_posteriors)
  # colnames(annotation_row) <- c("inferred clusters", "posterior")  
  # annotation_row$posterior <- as.factor(round(annotation_row$posterior,2))
  # annotation_row$`inferred clusters` <- as.factor(annotation_row$`inferred clusters`)
  
}

R <- dim(input_regions)[1] ## number of regions
M <- dim(input_CpG_data)[2] ## number of loci

#======================
# Plotting methylation data
#======================

## annotating the columns by region
reg_id <- unlist(sapply(1:R,function(x){rep(x,(input_regions[x,2]-input_regions[x,1])+1)}))
annotation_col <- as.matrix(reg_id,nrow=length(colnames(input_CpG_data)))
rownames(annotation_col) <- colnames(input_CpG_data)
annotation_col <- as.data.frame(annotation_col)
colnames(annotation_col) <- "regions"
annotation_col$regions <- as.factor(annotation_col$regions)

if (M <= 250) {
  print("plotting CpG based data")
  
  if(args$true_clusters == 1 && args$order_by_true == 1){
    tmp <- input_CpG_data[order(as.integer(true_cell_clusters$true_clusters)),]
  } else{
    tmp <- input_CpG_data[order(as.integer(inferred_cell_clusters$inferred_clusters)),]
  }
           
  pheatmap(tmp,cluster_rows = FALSE,cluster_cols=FALSE, annotation_row =annotation_row, cellwidth = 5,
           cellheight = 5,fontsize = 8, 
           main = paste0("CpG-based methylation data for ", args$name),
           gaps_row = index_gaps,fontsize_row=6,fontsize_col=4, 
           annotation_names_row = FALSE, annotation_names_col= FALSE,
           #gaps_col=(input_regions[,2][1:(R-1)] + 1),
           show_colnames=FALSE,
           annotation_col=annotation_col,
           filename = paste0(out_dir,"/",args$name,"_CpG_based_PLOT.pdf"))           
           
  rm(tmp)
  
  # ### attempt to plot hclust together
  # pheatmap(tmp,cluster_cols=FALSE, annotation_row = annotation_row,
  #          cellwidth = 5, cellheight = 5,
  #          fontsize = 8, main = "CpG-based methylation data",gaps_row = index_gaps,fontsize_row=6,fontsize_col=4, annotation_names_row = FALSE,
  #          show_colnames=FALSE,
  #          #gaps_col=(input_regions[,2][1:(R-1)]),
  #          annotation_col=annotation_col,
  #          filename = paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_CpG_based_PLOT.pdf"))
  # rm(tmp)
}

if (M > 250) {
  
  if (R == 1){
    
    print("plotting CpG based methylation data")
    
    if(args$true_clusters == 1 && args$order_by_true == 1){
      tmp <- input_CpG_data[order(as.integer(true_cell_clusters$true_clusters)),]
    } else{
      tmp <- input_CpG_data[order(as.integer(inferred_cell_clusters$inferred_clusters)),]
    }
             
    pheatmap(tmp,cluster_cols=FALSE, annotation_row = annotation_row,
             cluster_rows = FALSE,
             #cellwidth = 5, cellheight = 5,
             fontsize = 8, main = paste0("CpG based methylation data for ", args$name),
             gaps_row = index_gaps,fontsize_row=8,fontsize_col=4, 
             annotation_names_row = TRUE,annotation_names_col= FALSE,
             #gaps_col=(input_regions[,2][1:(R-1)] + 1),
             show_colnames=FALSE,
             annotation_col=annotation_col,
             filename = paste0(out_dir,"/",args$name,"_CpG_based_PLOT.pdf"))             
             
    rm(tmp)
    
  }
  
  if (R > 1){
    
    print("plotting CpG based methylation data")
    if(args$true_clusters == 1 && args$order_by_true == 1){
      tmp <- input_CpG_data[order(as.integer(true_cell_clusters$true_clusters)),]
    } else{
      tmp <- input_CpG_data[order(as.integer(inferred_cell_clusters$inferred_clusters)),]
    }
    
    ## in case we want to force some colors to be the same
    #colnames(annotation_row) <- "inferred_clusters"
    #annotation_row$inferred_clusters <- paste0("cl", annotation_row$inferred_clusters)
    #ann_colors = list(
    #inferred_clusters = c(cl0="blue", cl3="green",cl4="pink"))
    
    ## including the cell posterior probabilities in shades of grey
    #colnames(annotation_row) <- "inferred_clusters"
    #annotation_row$posterior <- paste0("post_", annotation_row$posterior) ### still need to make this more general
    #ann_colors = list(posterior = c(post_100="#111111", post_95="#333333",post_70="#666666", post_50="#777777", post_33="#888888", post_25="#999999" ))
             
    pheatmap(tmp,cluster_cols=FALSE, annotation_row = annotation_row,
             cluster_rows = FALSE,
             #cellwidth = 5, cellheight = 5,
             fontsize = 8, main = paste0("CpG based methylation data for ", args$name),
             gaps_row = index_gaps,fontsize_row=8,fontsize_col=4,
             annotation_names_row = FALSE,
             #annotation_colors = ann_colors,
             #annotation_names_col= TRUE,
             show_colnames=FALSE,
             #annotation_col=annotation_col,
             annotation_legend = TRUE,
             #legend_breaks = c(0,1),legend_labels = c("unmeth","meth"),
             filename = paste0(out_dir,"/",args$name,"_CpG_based_PLOT.pdf"))          
             
    rm(tmp)
    
    print("plotting region based mean methylation data")
    
    if(args$true_clusters == 1 && args$order_by_true == 1){
      tmp <- mean_meth_matrix[order(as.integer(true_cell_clusters$true_clusters)),]
    } else{
      tmp <- mean_meth_matrix[order(as.integer(inferred_cell_clusters$inferred_clusters)),]
    }
             
    pheatmap(tmp,cluster_cols=FALSE,
           cluster_rows=FALSE,
           annotation_row = annotation_row,
           #cellwidth = 5,cellheight = 5,
           fontsize = 8, main = paste0("Region-based mean methylation fraction data for ", args$name),
           gaps_row = index_gaps,fontsize_row=8,fontsize_col=6,
           annotation_names_row = FALSE,
           #annotation_colors = ann_colors,
           border_color=NA,show_colnames=FALSE,
           filename = paste0(out_dir,"/",args$name,"_region_based_PLOT.pdf"))             
             
    #filename = paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_region_based_PLOT.pdf"))
    rm(tmp)
    
    #=============================
    ### plotting CpG data from a subset of regions
    #=============================
    plot_subset <- 0
    if(plot_subset == 1){
      set.seed(14)
      input_regions <- input_regions[sort(sample(R,size=3)),] 
      input_CpG_data <- input_CpG_data[,c(input_regions[1,1]:input_regions[1,2],input_regions[2,1]:input_regions[2,2],input_regions[3,1]:input_regions[3,2])]
      
      R <- dim(input_regions)[1] ## number of regions
      M <- dim(input_CpG_data)[2] ## number of loci
      
      reg_id <- unlist(sapply(1:R,function(x){rep(x,(input_regions[x,2]-input_regions[x,1])+1)}))
      annotation_col <- as.matrix(reg_id,nrow=length(colnames(input_CpG_data)))
      rownames(annotation_col) <- colnames(input_CpG_data)
      annotation_col <- as.data.frame(annotation_col)
      colnames(annotation_col) <- "regions"
      annotation_col$regions <- as.factor(annotation_col$regions)
      
      
      if(args$true_clusters == 1 && args$order_by_true == 1){
        tmp <- input_CpG_data[order(as.integer(true_cell_clusters$true_clusters)),]
      } else{
        tmp <- input_CpG_data[order(as.integer(inferred_cell_clusters$inferred_clusters)),]
      }
               
      pheatmap(tmp,cluster_rows = FALSE,cluster_cols=FALSE, 
               #annotation_row =annotation_row, 
               cellwidth = 5, cellheight = 5,
               fontsize = 8, main = paste0("CpG-based methylation data for ", args$name),
               gaps_row = index_gaps,fontsize_row=6,fontsize_col=4, 
               annotation_names_row = FALSE,annotation_names_col= FALSE,
               show_colnames=FALSE,
               annotation_col=annotation_col,
               gaps_col =  (which(!duplicated(reg_id) == TRUE)[-1]-1),
               filename = paste0(out_dir,"/",args$name,"_subset_of_regions_CpG_based_PLOT.png"))
               #filename = paste0(sub(input_CpG_data_file,pattern=".tsv",replacement=""),"_subset_of_regions_CpG_based_PLOT.png"))               
      rm(tmp)
      
    }
    
    
    
  }
  
}

#======================
# Plotting CN data if available 
#======================

if (args$copy_number == 1 ){
  
  ## annotating the columns by chr, but it is currently not working very well
  annotation_col_chr <- as.matrix(colnames(input_CN_data),nrow=length(colnames(test)))
  rownames(annotation_col_chr) <- colnames(input_CN_data)
  annotation_col_chr <- as.data.frame(annotation_col_chr)
  colnames(annotation_col_chr) <- "chr"
  annotation_col_chr$chr <- paste0("chr",annotation_col_chr$chr)
  annotation_col_chr$chr <- as.factor(annotation_col_chr$chr)
  levels(annotation_col_chr$chr) <- paste0("chr",c(1:22,"X"))
  
  if (dim(input_CN_data)[2] < 250) {
    tmp <- input_CN_data[order(inferred_cell_clusters),]
             
    pheatmap(tmp,cluster_rows = FALSE,cluster_cols=FALSE, annotation_row = annotation_row, 
           #annotation_col = annotation_col_chr,
           cellwidth = 6,cellheight = 6,
           fontsize = 8, main = paste0("Copy number data for ", args$name),
           gaps_row = index_gaps,fontsize_row=6,fontsize_col=6, annotation_names_row = FALSE,
           filename = paste0(out_dir,"/",args$name,"_CN_PLOT.pdf"))             
    rm(tmp)
  }
  
  
  tmp <- input_CN_data[order(inferred_cell_clusters),]
           
  pheatmap(tmp,cluster_rows = FALSE,cluster_cols=FALSE, annotation_row = annotation_row,
             fontsize = 8, main = paste0("Copy number data for ", args$name),
             gaps_row = index_gaps,fontsize_row=6,fontsize_col=4,
             border_color=NA, annotation_names_row = FALSE,show_colnames=FALSE,
             annotation_col=annotation_col_chr,
             #annotation_colors = ann_colors,
             filename = paste0(out_dir,"/",args$name,"_noLines_CN_PLOT.pdf"))           
           
  rm(tmp)
  
  
}  


