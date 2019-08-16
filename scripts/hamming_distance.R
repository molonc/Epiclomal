
#======================
# arguments
#======================

suppressMessages(library(argparse))
library(pheatmap)
# create parser object
parser <- ArgumentParser()
parser$add_argument("--output_file", type="character", help="Path to the output file for the Bayesian estimates")
parser$add_argument("--true_epigenotype_file", type="character", help="Path to true epigenotypes")
parser$add_argument("--true_membership_file", type="character", help="Path to true cell membership")
parser$add_argument("--estimated_epigenotype_file", type="character",help="Path to estimated epigenotypes")
parser$add_argument("--estimated_membership_file", type="character", help="Path to estimated cell membership")
parser$add_argument("--methylation_file", type="character", help="Path to the methylation file")
parser$add_argument("--regions_file", type="character", help="Path to the region file")
args <- parser$parse_args()

# print(args)

outfile_est <- args$output_file

true_Z_file <- args$true_membership_file
true_epigenotypes_file <- args$true_epigenotype_file

estimated_Z_file <- args$estimated_membership_file
estimated_epigenotypes_file <- args$estimated_epigenotype_file
methylation_file <- args$methylation_file
regions_file <- args$regions_file

# to comment
# outfile_est <- "results_basic/DIC_LINE_ELBOW_gainthr0.05_0.02/all_hdist_bestrun_basic.tsv"
# estimated_Z_file <- "epiclomal_basic//20/cluster_MAP.tsv.gz"
# estimated_epigenotypes_file <- "epiclomal_basic//20/genotype_MAP.tsv.gz"
# true_Z_file <- "synthetic/data/true_clone_membership.tsv.gz"
# true_epigenotypes_file <- "synthetic/data/true_clone_epigenotypes.tsv.gz"
# methylation_file <- "synthetic/data/data_incomplete.tsv.gz"

# to comment
# outfile_est <- "results_region/DIC_LINE_ELBOW_gainthr0.05_0.02/all_hdist_bestrun_region.tsv"
# estimated_Z_file <- "epiclomal_region/20/cluster_MAP.tsv.gz"
# estimated_epigenotypes_file <- "epiclomal_region/20/genotype_MAP.tsv.gz"
# true_Z_file <- "synthetic/data/true_clone_membership.tsv.gz"
# true_epigenotypes_file <- "synthetic/data/true_clone_epigenotypes.tsv.gz"
# methylation_file <- "synthetic/data/data_incomplete.tsv.gz"
# regions_file <- "synthetic/data/regions_file.tsv.gz"

true_Z <- read.csv(true_Z_file,sep='\t')
true_epi <- read.csv(true_epigenotypes_file,sep='\t',header=TRUE)
# the number of rows is the number of true clusters

# to erase
#true_epi <- cbind(true_epi[,1], true_epi[,4924:5116])
# print("True epi")
# print(true_epi)

outfile_naive <- paste0(outfile_est, ".naive.tsv")
outfile_est_corr <- paste0(outfile_est, ".corr.tsv")

true_epi.tmp <- as.matrix(true_epi[,-1])
data_true <- (t(sapply(true_Z[,2],function(x){return(true_epi.tmp[which(true_epi[,1]==x),])}))) ### obtaining the true epigenotype for each cell

estimate_Z <- read.csv(estimated_Z_file,sep='\t')
estimate_epi <- read.csv(estimated_epigenotypes_file,sep='\t',header=TRUE)

# TO erase
#estimate_epi <- cbind(estimate_epi[,1],estimate_epi[,4924:5116])
# print("Estimate epi")
# print(estimate_epi)

#########
# compute the cellxCpG matrix of estimated CpGs
estimate_epi.tmp <- as.matrix(estimate_epi[,-1])
data_estimate <- (t(sapply(estimate_Z[,2],function(x){return(estimate_epi.tmp[which(estimate_epi[,1]==x),])})))

# this will be the corrected version of the data estimate
data_estimate_corr <- data_estimate
#########

meth_data <- read.csv(methylation_file,sep='\t',header=TRUE,row.names=1)
# to erase
#meth_data <- cbind(meth_data[,1],meth_data[,4924:5116])

regions <- read.csv(regions_file,sep='\t',header=TRUE,row.names=1)
num_regions <- dim(regions)[1]

incomplete_data <- meth_data

# impute the missing meth data the naive way
# we used to traverse by CpG, but now we traverse by region and then by CpG
#for (j in 1:dim(meth_data)[2]) {

# check if the region starts from 0 or 1
num_pred_clusters <- length(unique(estimate_Z[,2]))
add <- 0
if (regions[1,1] == 0) {
    add <- 1
}

cluster_set <- sort(unique(estimate_Z[,2]))

for (r in 1:num_regions) {
    rstart <- regions[r,1]+add
    rend <- regions[r,2]+add
    # check if this a variable region using Pearson correlation
    variable_region <- FALSE
    if (num_pred_clusters > 1) {
        for (cid1 in seq(1,length(cluster_set)-1)) {
            for (cid2 in seq(cid1+1,length(cluster_set))) {
                vector1 <- as.numeric(estimate_epi[estimate_epi$cluster_id==cluster_set[cid1],rstart:rend])
                vector2 <- as.numeric(estimate_epi[estimate_epi$cluster_id==cluster_set[cid2],rstart:rend])
                di <- sum(abs(vector1-vector2))/length(vector1)
                print(paste0("Region ", r, " cluster1 ", cluster_set[cid1], " cluster2 ", cluster_set[cid2], " di ", di))
                if (di > 0.5) {  # large distance
                    variable_region <- TRUE
                }
            }
        }
    }
    for (j in seq(rstart, rend)) {
        # traverse by cluster
        # print (paste0("Region ", r, " , CpG ", j))
        for (c in unique(estimate_Z[,2])) {
            # find all the cells in this cluster
            val <- floor(median(meth_data[estimate_Z[,2]==c,j],na.rm=TRUE))
            if (is.na(val)) {  # if none of the cells in this cluster had an observed value at this CpG
                # replace the naive matrix with the median of all cells
                val <- floor(median(meth_data[,j],na.rm=TRUE))
                if (is.na(val)) {
                    val <- 0
                }
                # also correct the estimated matrix only in this case and if the region is not variable
                # may have to check that this is a region that is mostly similar across clusters
                if (variable_region == FALSE) {
                    vec_corr <- meth_data[estimate_Z[,2]==c,j]
                    data_estimate_corr[estimate_Z[,2]==c,j] <- replace(vec_corr,is.na(vec_corr),val)
                }
            }
            vec <- meth_data[estimate_Z[,2]==c,j]
            meth_data[estimate_Z[,2]==c,j] <- replace(vec,is.na(vec),val)
        }
    }
}



############

compute_and_save_hd <- function (data_true, data_estimate, outfile) {
    ncells <- dim(data_true)[1]
    nloci <- dim(data_true)[2]
    hamming_distance_per_cell <- NULL

    for(i in 1:ncells){
      hd <- sum(data_true[i,] != data_estimate[i,])/nloci
      hamming_distance_per_cell <- c(hamming_distance_per_cell,hd)
    }

    # inferring the methylation profiles not from the G matrix, but from the clustering result
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
}


plot_cell_matrix <- function (matrix, file) {
    # get just a part of the matrix
    matrix <- matrix[1:40,4950:5100]
    pheatmap(matrix,cluster_rows = FALSE,cluster_cols=FALSE,
             annotation_names_row = FALSE, annotation_names_col= FALSE,
             cex = 0.8, cellwidth = 6, cellheight = 4,
             filename = file)
}

########


compute_and_save_hd (data_true, data_estimate, outfile_est)
compute_and_save_hd (data_true, data_estimate_corr, outfile_est_corr)
compute_and_save_hd (data_true, meth_data, outfile_naive)


#plot_cell_matrix (data_true, paste0(outfile_est, ".data_true_by_cell.pdf"))
#plot_cell_matrix (data_estimate, paste0(outfile_est, ".data_estimated_by_cell.pdf"))
#plot_cell_matrix (data_estimate_corr, paste0(outfile_est, ".data_estimated_corr_by_cell.pdf"))
#plot_cell_matrix (meth_data, paste0(outfile_est, ".data_naive_by_cell.pdf"))
#plot_cell_matrix (incomplete_data, paste0(outfile_est, ".data_original_by_cell.pdf"))

