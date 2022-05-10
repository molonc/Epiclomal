# Calculates overall stats in all cells and final regions for Epiclomal

#======================
# libraries
#======================
# .libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos5/3.2", "/clusterapp/software/linux-x86_64-centos5/R-3.2.3/lib64/R/library"))

# .libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.5", "/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.5"))

suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(plyr))
suppressMessages(library(pheatmap))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--data_ID", type="character", help="Some identification for the data, e.g., Smallwood2014_CGI")
parser$add_argument("--path_post_processed_CpG_data", type="character",help="Path to the directory containing the methylation data for each CpG in the reference genome for the regions considered. Each cell has its own tsv file")
parser$add_argument("--path_stats_region_data", type="character",help="Path to the directory containing the region based stats for each cell")

# Removing this argument, we'll do the 0 case when num_cells_cutoff is 0
# parser$add_argument("--type_cutoff", type="double", default=0, help="type of cutoff used to obtain filtered data (final set of regions), 0: by average missing proportion and 1: by number of cells with no more than a certain missing proportion")

parser$add_argument("--num_cells_cutoff", type="double", default=5, help="keeping regions with at least this number of cells with missing proportion < miss_prop_cutoff. If this number is 0, then it uses the average missing proportion for all cells")
parser$add_argument("--miss_prop_cutoff", type="double", default=0.95, help="missing proportion cutoff to obtain filtered data")
#parser$add_argument("--IQR_cutoff", type="double", default=0.0, help="keeping regions with IQR of region-based mean methylation across cells greater than the cutoff")
parser$add_argument("--nloci_cutoff", type="double", default=0.0, help="keeping regions with the highest IQR of region-based mean methylation across cells till number of loci is < nloci_cutoff. If nloci_cutoff = 0.0 it means this filtering is not applied")

## A parameter that combines the 3 parameters above into 1 argument separated by _
#parser$add_argument("--all_cutoffs", type="character", default=NULL, help="<num_cells_cutoff>_<miss_prop_cutoff>_<IQR_cutoff>. If given, it overwrites the other 3 arguments")

parser$add_argument("--all_cutoffs", type="character", default=NULL, help="<num_cells_cutoff>_<miss_prop_cutoff>_<nloci_cutoff>. If given, it overwrites the other 3 arguments")

parser$add_argument("--filter_regions_same_meth", type="character", default=NULL, help="if 1, not NULL, we will filter out regions with same methylation across cells")
parser$add_argument("--same_meth_cutoff", type="double", default=0.0, help="if 0.05 it means we allow 5% of regions (with data) to be different across cells, default is zero")


parser$add_argument("--plot_heatmap_unfiltered", type="double", default=1, help=" 1 for plotting and 0 for not plotting the heatmaps for unfiltered data")
parser$add_argument("--plot_heatmap_filtered", type="double", default=1, help=" 1 for plotting and 0 for not plotting the heatmaps for filtered data")

args <- parser$parse_args()

if (!is.null(args$all_cutoffs)) {
  tmp <- as.numeric(t(sapply(strsplit(args$all_cutoffs, split ="_"),function(x){return(x[1:3])})))
  #print(tmp)
  args$num_cells_cutoff = tmp[1]
  args$miss_prop_cutoff = tmp[2]
  #args$IQR_cutoff = tmp[3]
  args$nloci_cutoff = tmp[3]
}

if (!is.null(args$output_directory)) {
  dir.create(args$output_directory, showWarnings = FALSE, recursive=TRUE)
}


print(args)

#===========================
# auxiliary functions
#===========================

extraction_region_info_f <- function(region_info,type){  ### type should be "filtered" or "unfiltered"

  tmp <- region_info

  print(dim(tmp))

  print(head(tmp))

  tmp$region_id <- factor(tmp$region_id,levels=as.character(unique(tmp$region_id)))

  number_regions_single_CpG <- sum(tmp$region_cpgNum == 1)

  print(number_regions_single_CpG)

  pdf(file.path(outdir, paste0("hist_numCpGs_per_region_",type,"_",args$data_ID,".pdf")))
  hist(tmp$region_cpgNum,main="Number of CpGs per region",xlab="Number of CpGs")
  dev.off()

  pdf(file.path(outdir, paste0("boxplot_numCpGs_per_region_",type,"_",args$data_ID,".pdf")))
  boxplot(tmp$region_cpgNum,main="Number of CpGs per region")
  dev.off()

  pdf(file.path(outdir, paste0("hist_length_per_region_",type,"_",args$data_ID,".pdf")))
  hist(tmp$region_length,main="Region length - bp",xlab="bp")
  dev.off()

  pdf(file.path(outdir, paste0("boxplot_length_per_region_",type,"_",args$data_ID,".pdf")))
  boxplot(tmp$region_length,main="Region length - bp")
  dev.off()

  ### Saving as a table with some stats on number of CpGs and region length

  length_of_regions <- as.matrix(summary(tmp$region_length))

  print("Region of biggest size")
  print(tmp[tmp$region_length == max(tmp$region_length),])

  number_of_CpGs_per_region <- as.matrix(summary(tmp$region_cpgNum))

  stats <- as.matrix(c("min","1st_qu","median","mean","3rd_qu","max"))

  summary_info <- data.frame(cbind(stats,length_of_regions,number_of_CpGs_per_region))

  colnames(summary_info) <- c("stats","length_of_regions","number_of_CpGs_per_region")

  print(summary_info)

  write.table(summary_info, file.path(outdir, paste0("region_summary_info_",type,"_data",args$data_ID,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)

  rm(tmp)

  ### Some stats on distance between CpGs

  return(number_regions_single_CpG)

}

#######################################################

region_distance_CpGs_f <- function(CpG_based_data,type){ ### CpG_based_data is the region info from one of the cell files that contain also the CpG coordinates

  ### Some plots on distance between CpGs

  tmp <- CpG_based_data

  tmp$region_id <- factor(tmp$region_id,levels=as.character(unique(tmp$region_id)))

  mean_distance_CpG_per_region <- ddply(tmp, .(region_id), summarise,dist_mean = mean(diff(CpG_start),na.rm=TRUE))

  mean_distance_CpG_per_region$dist_mean[which(mean_distance_CpG_per_region$dist_mean == "NaN")] <- NA ### NAs correspond to regions with only one CpG

  pdf(file.path(outdir, paste0("hist_dist_CpGs_",type,"_",args$data_ID,".pdf")))
  hist(mean_distance_CpG_per_region$dist_mean,main="Mean distance between CpGs across regions",xlab="Mean distance")
  dev.off()

  pdf(file.path(outdir, paste0("boxplot_dist_CpGs_",type,"_",args$data_ID,".pdf")))
  boxplot(mean_distance_CpG_per_region$dist_mean,main="Mean distance between CpGs across regions")
  dev.off()

  x1 <- tmp$region_cpgNum[!duplicated(tmp$region_id)]
  x2 <- mean_distance_CpG_per_region$dist_mean
  bigdf = as.data.frame(cbind(x1,x2))

  p = ggplot(bigdf,aes(x=x1,y=x2))  +
    xlab("Number of CpG sites in regions") +
    ylab("Mean distance between CpGs across regions")
  p2 = p + geom_point(alpha = 0.01, colour="blue") +
    theme_bw()

  ggsave(file.path(outdir, paste0("scatter_dist_numCpGs_",type,"_",args$data_ID,".pdf")))

}

###################

num_cells_miss_function <- function(x,cutoff){

  y <- sum(x <= cutoff)

  return(y)

}

##############################################

same.meth.f <- function(x,same_cutoff){
  if((sum(!is.na(x))*same_cutoff) < 1){
    cutoff <- 0
  }else{
    cutoff <- ceiling(sum(!is.na(x))*same_cutoff)}

  unique_meth <- length(unique(round(x[!is.na(x)],5)))  ### I am rounding up to 5 decimals

  if(unique_meth > cutoff){
    if(unique_meth == 1){
      same_meth <- TRUE
    } else {
      same_meth <- FALSE
    }
  }

  if(unique_meth <= cutoff){
    same_meth <- TRUE
  }

  return(same_meth)
}



#===============================
# end of auxiliary functions
#===============================

outdir <- args$output_directory

all_CpG_cell_files <- list.files(path=args$path_post_processed_CpG_data, pattern="*.tsv.gz")

number_cells <- length(all_CpG_cell_files) ## this will be used to construct a table with the info about the data set

print(number_cells)

all_stats_cell_files <- list.files(args$path_stats_region_data, pattern="*.tsv")


################################################################################################################################
### Extracting region info for unfiltered data (the same for all cells, so I will use the data from just one of the files) #####
################################################################################################################################

print("Extracting region info")

tmp0 <- read.csv(file.path(args$path_post_processed_CpG_data, all_CpG_cell_files[1]),sep="\t",header=TRUE)

tmp0 <- tmp0[,1:8]

num_loci_unfiltered <- dim(tmp0)[1]

r_index <- !duplicated(tmp0$region_id) ### we need at least more than one region

print("total number of regions - unfiltered data")
total_number_regions_unfiltered <- sum(r_index)
print(total_number_regions_unfiltered)

tmp <- tmp0[r_index,]

table_unfiltered_region_info <- tmp

chrs_in_order <- as.character(table_unfiltered_region_info$chr[!duplicated(table_unfiltered_region_info$chr)])

print(chrs_in_order)

rm(r_index)

number_regions_single_CpG <- extraction_region_info_f(region_info=tmp,type="unfiltered") ## returns a number and produces plots and tables

regions_average_number_CpGs <-  read.table(file=file.path(outdir, paste0("region_summary_info_","unfiltered","_data",args$data_ID,".tsv")),sep="\t",header=TRUE)[4,3]

print(regions_average_number_CpGs)

region_distance_CpGs_f(CpG_based_data = tmp0,type="unfiltered") ## produces plots

print("End of extracting region info")

#### End of region stats ####

######################################################################
### Checking average missing proportion per cell - unfiltered data ###
######################################################################

### TO TEST

if(args$nloci_cutoff > 1) {

  print("checking missing proportion per cell")

  miss_prop_per_cell <- numeric(length(all_CpG_cell_files))

  # load cached results in case script did not complete last time it was run
  start <- 1
  cached_data <- file.path(outdir, 'missing_prop_per_cell.Rda.gz')
  if (file.exists(cached_data)) {
    load(cached_data)
    print("loading cached missing proportion per cell")
  }

  if (start != FALSE) {
    print(paste('starting from cell', start))
    for (c in start:length(all_CpG_cell_files)) {

      print(all_CpG_cell_files[c])

      tmp <- read.csv(file.path(args$path_post_processed_CpG_data, all_CpG_cell_files[c]),sep="\t",header=TRUE)

      if(c == 1){
        num_loci_unfiltered <- dim(tmp)[1]

        # CpG_nodata <- as.vector(matrix(0,num_loci_unfiltered,1))

        CpG_nodata_df <- as.data.frame(as.character(tmp$region_id))

        colnames(CpG_nodata_df) <- c("region_id")

        CpG_nodata_df$region_id <- as.character(CpG_nodata_df$region_id)

        CpG_nodata_df$CpG_nodata <- as.vector(matrix(0,num_loci_unfiltered,1))

      }

      CpG_nodata_df$CpG_nodata[!is.na(tmp$meth_frac)] <- 1

      #print(head(CpG_nodata_df$CpG_nodata[!is.na((tmp$meth_frac))]))

      print(paste("CpG sites with no data", sum(CpG_nodata_df$CpG_nodata == 0)))

      miss_prop_per_cell[c] <- sum(is.na(tmp$meth_frac))/num_loci_unfiltered

      start <- start + 1

      save(num_loci_unfiltered, CpG_nodata_df, miss_prop_per_cell, start, file = cached_data, compress = "gzip")
    }
    start <- FALSE
    save(num_loci_unfiltered, CpG_nodata_df, miss_prop_per_cell, start, file = cached_data, compress = "gzip")
  } else { print("data already processed, loaded from cache") }

  print("Number of loci with no data across all cells = ")
  print(sum(CpG_nodata_df$CpG_nodata == 0))

  ave_miss_prop <- mean(miss_prop_per_cell)

  pdf(file.path(outdir, paste0("boxplot_miss_prop_per_cell_unfiltered_",args$data_ID,".pdf")))
  boxplot(miss_prop_per_cell,main="Missing proportion per cell")
  dev.off()

  print("END of checking missing proportion per cell")

  ##########################################################
  ### Saving a table with some info for unfiltered data ####
  ##########################################################

  info_unfiltered <- as.matrix(c(num_loci_unfiltered,number_cells,total_number_regions_unfiltered,number_regions_single_CpG,ave_miss_prop))

  rownames(info_unfiltered) <- c("number of loci","number of cells","number of regions","number of regions containing only 1 CpG", "average missing proportion")

  print("Table with some info for unfiltered data")
  print(info_unfiltered)

  write.table(info_unfiltered, file.path(outdir, paste0("unfiltered_data_info_",args$data_ID,".tsv")), sep="\t", quote=FALSE, col.names=FALSE, append=TRUE)

}


#####################################################################################
## Creating matrices with region-based data obtained in Step 3, that is, putting  ###
## stats from all cells together                                                  ###
#####################################################################################

print("Creating matrices with region-based info - IQR, mean methylation, missing proportion")

cell_stats <- read.csv(file.path(args$path_stats_region_data ,all_stats_cell_files[1]),sep="\t",header=TRUE)
num_regions <- length(cell_stats$region_id)

region_mean_meth <- matrix(, nrow = num_regions, ncol = length(all_stats_cell_files))
region_miss_prop <- matrix(, nrow = num_regions, ncol = length(all_stats_cell_files))
region_IQR_meth <- matrix(, nrow = num_regions, ncol = length(all_stats_cell_files))
cell_ID <- character(length(all_stats_cell_files))

start <- 1
cached_data <- file.path(outdir, 'region_based_stats.Rda.gz')
if (file.exists(cached_data)){
  print("loading cached region based stats")
  load(cached_data)
}
if (start != FALSE) {
  print(paste('starting from cell', start))

  for (f in start:length(all_stats_cell_files)) {

    cell_stats <- read.csv(file.path(args$path_stats_region_data ,all_stats_cell_files[f]),sep="\t",header=TRUE)

    print(as.character(unique(cell_stats$cell_id)))
    cell_ID[f] <- as.character(unique(cell_stats$cell_id))

    cell_stats$region_id <- factor(cell_stats$region_id,levels=as.character(unique(cell_stats$region_id)))

    region_mean_meth[,f] <- cell_stats$region_mean
    region_miss_prop[,f] <- cell_stats$region_miss
    region_IQR_meth[,f] <- cell_stats$region_IQR

    start <- start + 1

    save(cell_ID, region_mean_meth, region_miss_prop, region_IQR_meth, start, file = cached_data, compress = "gzip")
  }

  start <- FALSE
  save(cell_ID, region_mean_meth, region_miss_prop, region_IQR_meth, start, file = cached_data, compress = "gzip")

} else { print("data already processed, loaded from cache") }

print(length(cell_ID))
print(length(cell_stats$region_id))
print(dim(region_mean_meth))
print(dim(region_miss_prop))
print(dim(region_IQR_meth))

colnames(region_mean_meth) <- cell_ID
rownames(region_mean_meth) <- as.character(cell_stats$region_id)

colnames(region_miss_prop) <- cell_ID
rownames(region_miss_prop) <- as.character(cell_stats$region_id)

colnames(region_IQR_meth) <- cell_ID
rownames(region_IQR_meth) <- as.character(cell_stats$region_id)

no_data_function <- function(x){
  if(sum(x == 1) == length(x)){
    a <- TRUE
  } else {
    a <- FALSE
  }
  return(a)
}

### regions with no data across all cells
regions_no_data_index <- apply(region_miss_prop,1,no_data_function)
print("Number of regions with no data across all cells:")
total_regions_no_data <- sum(regions_no_data_index == TRUE)
print(total_regions_no_data)

more_info_unfiltered <- as.matrix(c(total_regions_no_data))

rownames(more_info_unfiltered) <- c("number of regions with no data")

write.table(more_info_unfiltered, file.path(outdir, paste0("unfiltered_data_info_",args$data_ID,".tsv")),sep="\t",quote=FALSE,col.names=FALSE,append=TRUE)


#############################################################################################
### Region based heatmaps for missing proportion, methylation IQR and mean methylation  #####
#############################################################################################

if(args$plot_heatmap_unfiltered == 1){

  print("plotting heatmaps unfiltered")

  ### Plotting IQR, mean methylation and missing proportion for each cell and each region with data in at least one cell

  pheatmap(t(region_mean_meth[regions_no_data_index==FALSE,]),cluster_rows = FALSE,cluster_cols=FALSE,fontsize = 8,
           fontsize_row=6,fontsize_col=4,show_colnames = FALSE,
           filename = file.path(outdir, paste0("heatmap_region_mean_meth_unfiltered_data_",args$data_ID,".pdf")))

  pheatmap(t(region_miss_prop[regions_no_data_index==FALSE,]),cluster_rows = FALSE,cluster_cols=FALSE,fontsize = 8,
           fontsize_row=6,fontsize_col=4,show_colnames = FALSE,
           filename = file.path(outdir, paste0("heatmap_region_miss_prop_unfiltered_data_",args$data_ID,".pdf")))

  pheatmap(t(region_IQR_meth[regions_no_data_index==FALSE,]),cluster_rows = FALSE,cluster_cols=FALSE,fontsize = 8,
           fontsize_row=6,fontsize_col=4,show_colnames = FALSE,
           filename = file.path(outdir, paste0("heatmap_region_IQR_meth_unfiltered_data_",args$data_ID,".pdf")))

}

print("End of creating matrices")

###########################################################
## Finding IQR of region mean methylation across cells ####
###########################################################

print("Finding IQR of region mean methylation across cells")

if (args$filter_regions_same_meth == 1){

  print("applying new filter") ## Applying the new filter so that regions with same methylation across all cells are eliminated

  same_meth <-as.matrix(apply(region_mean_meth,1,same.meth.f,same_cutoff=args$same_meth_cutoff))

  rownames(same_meth) <- rownames(region_mean_meth)

  write.table(same_meth,file=file.path(outdir, paste0("regions_same_mean_meth_",args$data_ID,".tsv")),row.names=TRUE,col.names=FALSE,sep="\t",quote=FALSE)

}

IQR_meth <-as.matrix(apply(region_mean_meth,1,function(x){IQR(x,na.rm=TRUE)}))

rownames(IQR_meth) <- rownames(region_mean_meth)

write.table(IQR_meth,file=file.path(outdir, paste0("IQR_mean_meth_region_",args$data_ID,".tsv")),row.names=TRUE,col.names=FALSE,sep="\t",quote=FALSE)

#####################################################
## Finding average missing proportion per region ####
#####################################################

print("Finding average missing proportion per region")

mean_miss_prop <- as.matrix(apply(region_miss_prop,1,function(x){mean(x,na.rm=TRUE)}))

rownames(mean_miss_prop) <- rownames(region_miss_prop)

pdf(file.path(outdir, paste0("hist_region_mean_miss_prop_unfiltered_",args$data_ID,".pdf")))
hist(mean_miss_prop,main="Mean missing proportion per region across cells",xlab="Mean missing proportion")
dev.off()

pdf(file.path(outdir, paste0("boxplot_region_mean_miss_prop_unfiltered_",args$data_ID,".pdf")))
boxplot(mean_miss_prop,main="Mean missing proportion per region across cells")
dev.off()

write.table(mean_miss_prop,file=file.path(outdir, paste0("mean_miss_prop_region_",args$data_ID,".tsv")),row.names=TRUE,col.names=FALSE,sep="\t",quote=FALSE)

####################################################################################################
## Finding how many cells with missing proportion <= 95% (or whatever the cutoff is) per region ####
####################################################################################################

number_cells_miss_cutoff <- as.matrix((apply(region_miss_prop,1,num_cells_miss_function,cutoff=args$miss_prop_cutoff)))

pdf(file.path(outdir, paste0("hist_region_number_cells_miss_cutoff_unfiltered_",args$data_ID,".pdf")))
hist(number_cells_miss_cutoff,main="# of cells with missing proportion < cutoff per region",xlab="Number of cells")
dev.off()

### Creating data.frame with region_id + miss_prop

tmp1 <- as.character(rownames(mean_miss_prop))
tmp2 <- as.vector(mean_miss_prop)

tmp3 <- data.frame(region_id = as.character(tmp1), miss_prop = tmp2)

mean_miss_prop <- tmp3

rm(tmp1,tmp2,tmp3)

### Creating data.frame with region_id + number_cells_miss_cutoff

tmp1 <- as.character(rownames(number_cells_miss_cutoff))
tmp2 <- as.vector(number_cells_miss_cutoff)

tmp3 <- data.frame(region_id = as.character(tmp1), number_cells_cutoff = tmp2)

number_cells_miss_cutoff <- tmp3

rm(tmp1,tmp2,tmp3)

### Checking number of regions of absolutely no data

print("Number of regions with no data")
regions_no_data <- sum((mean_miss_prop[,2] == 1))
print(regions_no_data)

print("Number of regions with zero cells meeting the miss_prop cutoff")
print(sum(number_cells_miss_cutoff[,2] == 0))

### Creating data frame for IQR_meth as well

tmp1 <- as.character(rownames(IQR_meth))
tmp2 <- as.vector(IQR_meth)

tmp3 <- data.frame(tmp1,tmp2)
colnames(tmp3) <- c("region_id","IQR")

tmp3$region_id <- as.character(tmp3$region_id)

IQR_meth <- tmp3

rm(tmp1,tmp2,tmp3)

# ### just doing some test to confirm that regions with miss_prop = 1 corresponds to regions where IQR = NA
#
# print("Tests:")
# print(sum((is.na(IQR_meth$IQR))))
#
# if(sum(is.na(IQR_meth$IQR)) == 0){
#
#   print(sum( IQR_meth$region_id != mean_miss_prop$region_id ))
#
# }
#
# if(sum(is.na(IQR_meth$IQR)) > 0){
#
#   print(sum( IQR_meth$region_id[is.na(IQR_meth$IQR)] == mean_miss_prop$region_id[(mean_miss_prop$miss_prop == 1)] ))
#
# }
#
# print("End of tests")


##########################
## Applying filters  #####
##########################

if (args$nloci_cutoff > 1 ){

  print("FILTERING BY NUMBER OF LOCI")

  regions_passed_miss_cutoff <- table_unfiltered_region_info
  dim(regions_passed_miss_cutoff)
  dim(IQR_meth)

  CpG_nodata_df$region_id <- factor(CpG_nodata_df$region_id,levels=as.character(regions_passed_miss_cutoff$region_id))

  new_number_CpGs <- ddply(CpG_nodata_df, .(region_id), summarise,number_CpGs = sum(CpG_nodata == 1)) ### 1 means data present for that CpG

  regions_passed_miss_cutoff$region_cpgNum <- new_number_CpGs$number_CpGs

  print(dim(regions_passed_miss_cutoff))
  print(dim(CpG_nodata_df))

  if ((args$filter_regions_same_meth == 1)) {

    print("Filtering regions with same methylation")

    index_same <- same_meth

    IQR_meth <- IQR_meth[!index_same,]

    mean_miss_prop <- mean_miss_prop[!index_same,]

    number_cells_miss_cutoff <- number_cells_miss_cutoff[!index_same,]

    region_mean_meth <- region_mean_meth[!index_same,]
    region_IQR_meth <- region_IQR_meth[!index_same,]
    region_miss_prop <- region_miss_prop[!index_same,]

    regions_passed_miss_cutoff <- regions_passed_miss_cutoff[!index_same,]

    print("End of applying new filter")

  }

  if(args$num_cells_cutoff == 0) {

    ###################################################################################
    ### finding regions with less than a certain amount of data missing across all  ###
    ### cells and then regions with certain IQR up to a certain number of loci      ###
    ###################################################################################

    print("Finding final set of regions by average missing proportion type of cutoff")

    index_miss <- which(mean_miss_prop$miss_prop <= args$miss_prop_cutoff)

    regions_passed_miss_cutoff <- regions_passed_miss_cutoff[index_miss,]

    number_CpGs_passed_miss_cutoff <- sum(regions_passed_miss_cutoff$region_cpgNum)

    #print(dim(regions_passed_miss_cutoff))

    print("Number of CpGs that passed miss_prop cutoff")

    print(number_CpGs_passed_miss_cutoff)

    if ((number_CpGs_passed_miss_cutoff + regions_average_number_CpGs) <= args$nloci_cutoff){

      IQR_cutoff_implemented <- 0

      print("Not enough loci left to apply IQR cut off ")

      IQR_meth_tmp <- IQR_meth[index_miss,]

      regions_passed_IQR_cutoff <- regions_passed_miss_cutoff

      number_CpGs_passed_IQR_cutoff <- sum(regions_passed_IQR_cutoff$region_cpgNum)

      regions_IQR_greater <- min(IQR_meth_tmp$IQR[length(IQR_meth_tmp$IQR)])

      save_table <- as.matrix(c(dim(regions_passed_IQR_cutoff)[1],number_CpGs_passed_IQR_cutoff,regions_IQR_greater,IQR_cutoff_implemented))

      rownames(save_table) <- c("number of regions that passed IQR/nloci cutoff","number of final CpGs","smallest region mean IQR","IQR cutoff implemented")

      print(save_table)

      write.table(save_table,file=file.path(outdir, paste0("extra_info_after_applying_cutoffs_",args$data_ID,".tsv")),row.names=TRUE,col.names=FALSE,sep="\t",quote=FALSE)

      write.table(IQR_meth_tmp$region_id,file=file.path(outdir, paste0("final_regions_",args$data_ID,".tsv")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

      index_IQR <- 1:dim(IQR_meth_tmp)[1]

    }

    if ((number_CpGs_passed_miss_cutoff + regions_average_number_CpGs)  > args$nloci_cutoff) {

      IQR_cutoff_implemented <- 1

      print("IQR CUTOFF IMPLEMENTED")

      IQR_meth_tmp <- IQR_meth[index_miss,]

      index_IQR_tmp <- order(IQR_meth_tmp$IQR,decreasing=TRUE)[1:(dim(IQR_meth_tmp)[1])]

      regions_passed_miss_cutoff <- regions_passed_miss_cutoff[index_IQR_tmp,] ### ordering by IQR

      IQR_meth_tmp <- IQR_meth_tmp[index_IQR_tmp,] ### ordering by IQR

      ### After ordering the regions by IQR keeping the first given number of loci (argument == nloci_cutoff)
      index_IQR <- (which(cumsum(regions_passed_miss_cutoff$region_cpgNum) <= args$nloci_cutoff))

      number_regions_cutoff <- length(index_IQR)

      regions_passed_IQR_cutoff <- regions_passed_miss_cutoff[index_IQR,]

      number_CpGs_passed_IQR_cutoff <- sum(regions_passed_IQR_cutoff$region_cpgNum)

      regions_IQR_greater <- IQR_meth_tmp[index_IQR,]$IQR[length(index_IQR)]

      save_table <- as.matrix(c(dim(regions_passed_IQR_cutoff)[1],number_CpGs_passed_IQR_cutoff,regions_IQR_greater,IQR_cutoff_implemented))

      rownames(save_table) <- c("number of regions that passed IQR/nloci cutoff","number of final CpGs","smallest region mean IQR","IQR cutoff implemented")

      print(save_table)

      write.table(save_table,file=file.path(outdir, paste0("extra_info_after_applying_cutoffs_",args$data_ID,".tsv")),row.names=TRUE,col.names=FALSE,sep="\t",quote=FALSE)

      final_regions_unordered <- IQR_meth_tmp$region_id[index_IQR]

      prev_IQR_meth_tmp <- IQR_meth[index_miss,]

      ### new index_IQR that will give ordered regions:
      index_IQR <- sort(as.numeric(sapply(final_regions_unordered,function(x){which(prev_IQR_meth_tmp$region_id==x)})))

      FINAL_reg <- prev_IQR_meth_tmp$region_id[index_IQR]

      print(head(FINAL_reg))

      print(length(FINAL_reg))

      write.table(FINAL_reg,file=file.path(outdir, paste0("final_regions_",args$data_ID,".tsv")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

    }
  } else {

    ####################################################################################################################
    ### finding regions with at least 5 (or other value) cells with less than 95% (or other value) of missing data   ###
    ### and then regions with certain IQR                                                                            ###
    ####################################################################################################################

    print("Finding final set of regions by using as cutoff number of cells with certain missing proportion")

    #print(dim(number_cells_miss_cutoff))

    #print(dim(regions_passed_miss_cutoff))

    index_miss <- which(number_cells_miss_cutoff$number_cells_cutoff >= args$num_cells_cutoff)

    regions_passed_miss_cutoff <- regions_passed_miss_cutoff[index_miss,]

    number_CpGs_passed_miss_cutoff <- sum(regions_passed_miss_cutoff$region_cpgNum)

    #print(number_CpGs_passed_miss_cutoff)
    #print(regions_average_number_CpGs)

    if ((number_CpGs_passed_miss_cutoff + regions_average_number_CpGs) <= args$nloci_cutoff) {

      IQR_cutoff_implemented <- 0

      print("Not enough loci left to apply IQR cut off ")

      IQR_meth_tmp <- IQR_meth[index_miss,]

      regions_passed_IQR_cutoff <- regions_passed_miss_cutoff

      number_CpGs_passed_IQR_cutoff <- sum(regions_passed_IQR_cutoff$region_cpgNum)

      regions_IQR_greater <- min(IQR_meth_tmp$IQR[length(IQR_meth_tmp$IQR)])

      save_table <- as.matrix(c(dim(regions_passed_IQR_cutoff)[1],number_CpGs_passed_IQR_cutoff,regions_IQR_greater,IQR_cutoff_implemented))

      rownames(save_table) <- c("number of regions that passed IQR/nloci cutoff","number of final CpGs","smallest region mean IQR","IQR cutoff implemented")

      print(save_table)

      write.table(save_table,file=file.path(outdir, paste0("extra_info_after_applying_cutoffs_",args$data_ID,".tsv")),row.names=TRUE,col.names=FALSE,sep="\t",quote=FALSE)

      write.table(IQR_meth_tmp$region_id,file=file.path(outdir, paste0("final_regions_",args$data_ID,".tsv")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

      index_IQR <- 1:dim(IQR_meth_tmp)[1]
    }

    if ((number_CpGs_passed_miss_cutoff + regions_average_number_CpGs) > args$nloci_cutoff){

      IQR_cutoff_implemented <- 1

      print("IQR CUTOFF IMPLEMENTED")

      IQR_meth_tmp <- IQR_meth[index_miss,]

      index_IQR_tmp <- order(IQR_meth_tmp$IQR,decreasing=TRUE)[1:(dim(IQR_meth_tmp)[1])]

      regions_passed_miss_cutoff <- regions_passed_miss_cutoff[index_IQR_tmp,] ### ordering by IQR

      IQR_meth_tmp <- IQR_meth_tmp[index_IQR_tmp,] ### ordering by IQR

      ### After ordering the regions by IQR keeping the first given number of loci (argument == nloci_cutoff)
      index_IQR <- which(cumsum(regions_passed_miss_cutoff$region_cpgNum) <= args$nloci_cutoff)

      number_regions_cutoff <- length(index_IQR)

      regions_passed_IQR_cutoff <- regions_passed_miss_cutoff[index_IQR,]

      number_CpGs_passed_IQR_cutoff <- sum(regions_passed_IQR_cutoff$region_cpgNum)

      regions_IQR_greater <- IQR_meth_tmp[index_IQR,]$IQR[length(index_IQR)]

      save_table <- as.matrix(c(dim(regions_passed_IQR_cutoff)[1],number_CpGs_passed_IQR_cutoff,regions_IQR_greater,IQR_cutoff_implemented))

      rownames(save_table) <- c("number of regions that passed IQR/nloci cutoff","number of final CpGs","smallest region mean IQR","IQR cutoff implemented")

      print(save_table)

      write.table(save_table,file=file.path(outdir, paste0("extra_info_after_applying_cutoffs_",args$data_ID,".tsv")),row.names=TRUE,col.names=FALSE,sep="\t",quote=FALSE)

      ### Ordering the regions by chr and start

      print("Ordering the regions by chr and start")

      ### Using 2nd way just as above for mean missing proportion

      final_regions_unordered <- IQR_meth_tmp$region_id[index_IQR]

      prev_IQR_meth_tmp <- IQR_meth[index_miss,]

      ### new index_IQR that will give ordered regions:
      index_IQR <- sort(as.numeric(sapply(final_regions_unordered,function(x){which(prev_IQR_meth_tmp$region_id==x)})))

      FINAL_reg <- prev_IQR_meth_tmp$region_id[index_IQR]

      print(head(FINAL_reg))

      write.table(FINAL_reg,file=file.path(outdir, paste0("final_regions_",args$data_ID,".tsv")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    }
  }

  #############################################################################################################
  ### saving the matrix with the region-based methylation for the final regions after applying the cutoffs ####
  #############################################################################################################

  print("saving the matrix with region-based methylation for final regions after applying cutoffs")

  new_region_mean_meth <- region_mean_meth[index_miss,]

  ##new_region_mean_meth <- new_region_mean_meth[index_IQR_tmp,]

  new_region_mean_meth <- new_region_mean_meth[index_IQR,]

  write.table(new_region_mean_meth,file=file.path(outdir, paste0("final_mean_meth_region_",args$data_ID,".tsv")),row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

  ################################################################################
  ### Plot heatmaps for filtered IQR, mean methylation and missing proportion ####
  ################################################################################

  new_region_IQR_meth <- region_IQR_meth[index_miss,]

  ##new_region_IQR_meth <- new_region_IQR_meth[index_IQR_tmp,]

  new_region_IQR_meth <- new_region_IQR_meth[index_IQR,]

  ### Filtered missing proportion matrix

  new_region_miss_prop <- region_miss_prop[index_miss,]

  ##new_region_miss_prop  <- new_region_miss_prop[index_IQR_tmp,]

  new_region_miss_prop  <- new_region_miss_prop[index_IQR,]

  ### Plots

  if(args$plot_heatmap_filtered == 1){

    ### Plotting IQR, mean methylation and missing proportion for each cell and each region for filtered data

    pheatmap(t(new_region_mean_meth),cluster_rows = FALSE,cluster_cols=FALSE,fontsize = 8,
             fontsize_row=6,fontsize_col=4,show_colnames = FALSE,
             filename = file.path(outdir, paste0("heatmap_region_mean_meth_filtered_data_",args$data_ID,".pdf")))

    pheatmap(t(new_region_miss_prop),cluster_rows = FALSE,cluster_cols=FALSE,fontsize = 8,
             fontsize_row=6,fontsize_col=4,show_colnames = FALSE,
             filename = file.path(outdir, paste0("heatmap_region_miss_prop_filtered_data_",args$data_ID,".pdf")))

    pheatmap(t(new_region_IQR_meth),cluster_rows = FALSE,cluster_cols=FALSE,fontsize = 8,
             fontsize_row=6,fontsize_col=4,show_colnames = FALSE,
             filename = file.path(outdir, paste0("heatmap_region_IQR_meth_filtered_data_",args$data_ID,".pdf")))

  }


  #########################################################
  ### Region based stats and plots for filtered data ######
  #########################################################

  filtered_regions <- IQR_meth_tmp[index_IQR,]$region_id

  print("Extracting region info for filtered data")

  tmp <- table_unfiltered_region_info

  tmp2 <- tmp[index_miss,]

  ##tmp2 <- tmp2[index_IQR_tmp,]

  tmp3 <- tmp2[index_IQR,]

  number_regions_single_CpG <- extraction_region_info_f(region_info=tmp3,type="filtered") ## returns a number and produces plots and tables
  print("Number of regions with only one CpG - filtered data")
  print(number_regions_single_CpG)
}


########################################################################
###### Another way: considering IQR values not number of loci     ######
########################################################################

if (args$nloci_cutoff <= 1 ) {

  print("FILTERING BY IQR")

  if (args$filter_regions_same_meth == 1) {

    print("Applying new filter")

    index_same <- same_meth

    IQR_meth <- IQR_meth[!index_same,]

    mean_miss_prop <- mean_miss_prop[!index_same,]

    number_cells_miss_cutoff <- number_cells_miss_cutoff[!index_same,]

    region_mean_meth <- region_mean_meth[!index_same,]

    region_IQR_meth <- region_IQR_meth[!index_same,]

    region_miss_prop <- region_miss_prop[!index_same,]

    print("End of applying new filter")

  }

  print("applying missing proportion and IQR filters")

  if(args$num_cells_cutoff == 0) {

    ###################################################################################
    ### finding regions with less than a certain amount of data missing across all  ###
    ### cells and then regions with certain IQR                                     ###
    ###################################################################################

    print("Finding final set of regions by average missing proportion type of cutoff")
    print(paste0("miss prop cutoff = ",args$miss_prop_cutoff))

    index_miss <- which(mean_miss_prop$miss_prop <= args$miss_prop_cutoff)

    print(dim(IQR_meth))
    IQR_meth_tmp <- IQR_meth[index_miss,]
    print(dim(IQR_meth_tmp))

    ### now for the remaining regions doing the IQR cutoff as well

    print(paste0("IQR cutoff = ",args$nloci_cutoff))
    index_IQR <- which(IQR_meth_tmp$IQR >= args$nloci_cutoff)

    write.table(IQR_meth_tmp[index_IQR,]$region_id,file=file.path(outdir, paste0("final_regions_",args$data_ID,".tsv")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)

  } else {
    ####################################################################################################################
    ### finding regions with at least 5 (or other value) cells with less than 95% (or other value) of missing data   ###
    ### and then regions with certain IQR                                                                            ###
    ####################################################################################################################

    print("Finding final set of regions by using as cutoff number of cells with certain missing proportion")

    ### First finding the set of regions with missing proportion smaller than cutoff

    index_miss <- which(number_cells_miss_cutoff$number_cells_cutoff >= args$num_cells_cutoff)
    IQR_meth_tmp <- IQR_meth[index_miss,]

    ### now for the remaining regions doing the IQR cutoff as well

    index_IQR <- which(IQR_meth_tmp$IQR >= args$nloci_cutoff)
    print(dim(IQR_meth_tmp[index_IQR,]))

    write.table(IQR_meth_tmp[index_IQR,]$region_id,file=file.path(paste0("final_regions_",args$data_ID,".tsv")),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
  }

  print("End of applying all cutoffs")

  ###########################################################################################################
  ### saving the matrix with the region-based methylation for the final regions after applying the cutoffs ##
  ###########################################################################################################
  print("saving the matrix with the region-based methylation for the final regions after applying cutoffs")

  print(dim(region_mean_meth))
  new_region_mean_meth <- region_mean_meth[index_miss,]

  new_region_mean_meth <- new_region_mean_meth[index_IQR,]
  print(dim(new_region_mean_meth))

  write.table(new_region_mean_meth,file=file.path(outdir, paste0("final_mean_meth_region_",args$data_ID,".tsv")),row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

  ################################################################################
  ### Plot heatmaps for filtered IQR, mean methylation and missing proportion ####
  ################################################################################

  ### Filtered IQR matrix

  print(dim(region_IQR_meth))
  new_region_IQR_meth <- region_IQR_meth[index_miss,]
  new_region_IQR_meth <- new_region_IQR_meth[index_IQR,]
  print(dim(new_region_IQR_meth))

  ### Filtered missing proportion matrix

  print(dim(region_miss_prop))
  new_region_miss_prop <- region_miss_prop[index_miss,]
  new_region_miss_prop  <- new_region_miss_prop[index_IQR,]
  print(dim(new_region_miss_prop))

  ### Plots

  if(args$plot_heatmap_filtered == 1){

    print("plotting heatmaps after filtering cutoffs")

    ### Plotting IQR, mean methylation and missing proportion for each cell and each region for filtered data

    pheatmap(t(new_region_mean_meth),cluster_rows = FALSE,cluster_cols=FALSE,fontsize = 8,
             fontsize_row=6,fontsize_col=4,show_colnames = FALSE,
             filename = file.path(outdir, paste0("heatmap_region_mean_meth_filtered_data_",args$data_ID,".pdf")))

    pheatmap(t(new_region_miss_prop),cluster_rows = FALSE,cluster_cols=FALSE,fontsize = 8,
             fontsize_row=6,fontsize_col=4,show_colnames = FALSE,
             filename = file.path(outdir, paste0("heatmap_region_miss_prop_filtered_data_",args$data_ID,".pdf")))

    pheatmap(t(new_region_IQR_meth),cluster_rows = FALSE,cluster_cols=FALSE,fontsize = 8,
             fontsize_row=6,fontsize_col=4,show_colnames = FALSE,
             filename = file.path(outdir, paste0("heatmap_region_IQR_meth_filtered_data_",args$data_ID,".pdf")))

  }
}
print("DONE")
