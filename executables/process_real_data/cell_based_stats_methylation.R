#======================
# libraries
#======================
# .libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos5/3.2", "/clusterapp/software/linux-x86_64-centos5/R-3.2.3/lib64/R/library"))

suppressMessages(library(argparse))
suppressMessages(library(plyr))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--data_ID", type="character", help="Some identification for the data, e.g., Smallwood2014_CGI")
parser$add_argument("--post_processed_CpG_data_file", type="character",help="Path to cell file containing the methylation data for each CpG in the reference genome for the regions considered")
parser$add_argument("--cell_ID", type="character", help="Cell ID")

args <- parser$parse_args()

print(args)

# ============================
# Auxiliary Functions
# ============================

miss.prop.function <- function(x){

 prop <-  sum(is.na(x))/length(x)

 return(prop)
}


ave.cov.function <- function(x1,x2){

  ave_cov <-  mean((x1+x2),na.rm=TRUE)

  return(ave_cov)
}

# ============================
# End of Auxiliary Functions
# ============================

outdir <- args$output_directory

# Finding the input file from the path and cell_id
input_file <- args$post_processed_CpG_data_file
print("Input file for this cell is")
print(input_file)

tmp <- read.csv(input_file,sep="\t",header=TRUE)

print(dim(tmp))
  ### using a smaller set when testing
  # tmp <- tmp[1:1000,]

tmp$region_id <- factor(tmp$region_id,levels=as.character(unique(tmp$region_id)))

##############################################################
## Calculating various statistics for each region in a cell ##
##############################################################

stat_tmp <- ddply(tmp, .(region_id), summarise,region_mean = mean(meth_frac,na.rm=TRUE),region_median = median(meth_frac,na.rm=TRUE),region_IQR = IQR(meth_frac,na.rm=TRUE),
                  region_miss = miss.prop.function(meth_frac),region_cov = ave.cov.function(count_meth,count_unmeth))

stat_tmp$region_mean[which(stat_tmp$region_mean == "NaN")] <- NA
stat_tmp$region_IQR[which(stat_tmp$region_IQR == "NaN")] <- NA
stat_tmp$region_miss[which(stat_tmp$region_miss == "NaN")] <- NA
stat_tmp$region_cov[which(stat_tmp$region_cov == "NaN")] <- NA
stat_tmp$region_median[which(stat_tmp$region_median == "NaN")] <- NA

cell_id <- as.character(tmp$cell_id[!duplicated(tmp$region_id)])

stat_tmp <- cbind(stat_tmp,cell_id)

print(head(stat_tmp))

#####################################
### Saving the region based stats ###
#####################################

write.table(stat_tmp,file=file.path(outdir, paste0("stats_region_",args$data_ID,"_",args$cell_ID,".tsv")), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

print("DONE")

