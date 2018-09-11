
#======================
# arguments
#======================

suppressMessages(library(argparse))
# create parser object
parser <- ArgumentParser()
parser$add_argument("--output_file", type="character", help="Path to the output file") 
parser$add_argument("--true_epigenotype_file", type="character", help="Path to true epigenotypes") 
parser$add_argument("--true_membership_file", type="character", help="Path to true cell membership") 
parser$add_argument("--estimated_epigenotype_file", type="character",help="Path to estimated epigenotypes")
parser$add_argument("--estimated_membership_file", type="character", help="Path to estimated cell membership") 
args <- parser$parse_args() 

# print(args)

outfile <- args$output_file

true_Z_file <- args$true_membership_file
true_epigenotypes_file <- args$true_epigenotype_file

estimated_Z_file <- args$estimated_membership_file
estimated_epigenotypes_file <- args$estimated_epigenotype_file

#estimated_Z_file <- "~/Documents/hamming_distance/true_clone_membership.tsv"
#estimated_epigenotypes_file <- "~/Documents/hamming_distance/true_clone_epigenotypes.tsv"
#true_Z_file <- "~/Documents/hamming_distance/true_clone_membership.tsv"
#true_epigenotypes_file <- "~/Documents/hamming_distance/true_clone_epigenotypes.tsv"

true_Z <- read.csv(true_Z_file,sep='\t')
true_epi <- read.csv(true_epigenotypes_file,sep='\t',header=TRUE)

true_epi.tmp <- as.matrix(true_epi[,-1])
data_true <- (t(sapply(true_Z[,2],function(x){return(true_epi.tmp[which(true_epi[,1]==x),])}))) ### obtaining the true epigenotype for each cell

estimate_Z <- read.csv(estimated_Z_file,sep='\t')
estimate_epi <- read.csv(estimated_epigenotypes_file,sep='\t',header=TRUE)

estimate_epi.tmp <- as.matrix(estimate_epi[,-1])
data_estimate <- (t(sapply(estimate_Z[,2],function(x){return(estimate_epi.tmp[which(estimate_epi[,1]==x),])})))

ncells <- dim(data_true)[1]
nloci <- dim(data_true)[2]
hamming_distance_per_cell <- NULL

for(i in 1:ncells){
  hd <- sum(data_true[i,] != data_estimate[i,])/nloci
  hamming_distance_per_cell <- c(hamming_distance_per_cell,hd)
  
}

print("number of cells")
print(ncells)
print("number of CpGs")
print(nloci)

hd_stats <- t(as.matrix(summary(hamming_distance_per_cell)))
iqr <- IQR(hamming_distance_per_cell)
hd_stats <- cbind(hd_stats,iqr)
colnames(hd_stats) <- c("min","1stQu","median","mean","3rdQu","max","IQR")

print("Stats for the relative hamming distances across cells")
print(hd_stats)

### if we just want the median and IQR
#hd_stats <-t(as.matrix(c(median(hamming_distance_per_cell),IQR(hamming_distance_per_cell))))
#print("median cell relative hamming distance and IQR")
#colnames(hd_stats) <- c("median","IQR")


write.table (hd_stats, file=outfile, sep="\t", row.names=FALSE)


