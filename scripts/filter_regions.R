suppressMessages(library(argparse))
suppressMessages(library(pheatmap))

parser <- ArgumentParser()

parser$add_argument("--mean_methylation_file", default="None", type="character", help="file path contianing mean methylation per region")
parser$add_argument("--true_file", default=NULL, type="character", help="Path to true cluster information")
parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--coef_threshold", default=0.95, type="double",help="Threshold of Pearson correlation coefficient")
parser$add_argument("--mean_diff_threshold", default=0.05, type="double",help="Threshold of max difference between mean region cluster methylation")


args <- parser$parse_args()

print(args)

mean_meth_file <- args$mean_methylation_file
true_file <- args$true_file
outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)

coef_t <- args$coef_threshold

load_mean_meth <- function(filename) {
    mean_meth_matrix <- read.csv(filename, sep = "\t", header = TRUE)
    mean_meth_matrix <- t(mean_meth_matrix)
    return(mean_meth_matrix)
}

set_annotation_row <- function(data) {
    annotation_row <- as.data.frame(data$epigenotype_id)
    colnames(annotation_row) <- "true_cluster"
    annotation_row$true_cluster <- as.factor(annotation_row$true_cluster)
    rownames(annotation_row) <- data$cell_id
    return(annotation_row)
}

set_index_gaps <- function(data) {
    index <- 1:dim(data)[1]
    index_gaps <- index[!duplicated(data[order(data),])] - 1
    index_gaps <- index_gaps[which(index_gaps != 0)]
    return(index_gaps)
}

find_correlated_regions <- function(pearson_corr, coef_t) {
    correlated_regions <- NULL
    for (r in 1:(nrow(pearson_corr) - 1)) {
        if (rownames(pearson_corr)[r] %in% correlated_regions) {
            next
        }
        for (c in (r+1):nrow(pearson_corr)) {
            if (colnames(pearson_corr)[c] %in% correlated_regions) {
                next
            }
            if (!is.na(pearson_corr[r,c]) && ((pearson_corr[r,c] > coef_t) || (pearson_corr[r,c] < (-1 * coef_t)))) {
                correlated_regions <- append(correlated_regions, colnames(pearson_corr)[c])
            }
        }
    }
    return(correlated_regions)
}

find_low_v_regions <- function(data, true_clusters) {
    data <- t(data)
    epi_means <- NULL
    for (epi in unique(true_clusters$epigenotype_id)) {
        epi_means <- cbind(epi_means, rowMeans(data[,true_clusters$epigenotype_id == epi], na.rm = TRUE))
    }
    epi_means <- as.data.frame(epi_means)
    epi_means$max_diff <- apply(epi_means, 1, function(x) diff(range(x, na.rm = TRUE)))

    low_variance_regions <- rownames(epi_means[epi_means$max_diff < 0.05,])

    return(low_variance_regions)
}

main <- function(mean_meth_file, true_file, output_directory, coef_threshold) {
    mean_meth_matrix <- load_mean_meth(mean_meth_file)
    print(paste("Number of cells", dim(mean_meth_matrix)[1]))
    print(paste("Number of regions", dim(mean_meth_matrix)[2]))

    print("Calculating Pearson correlation coefficient between regions")
    pearson_corr <- cor(x=mean_meth_matrix, method="pearson", use="pairwise.complete.obs")
    print("Finding and removing highly correlated regions")
    correlated_regions <- find_correlated_regions(pearson_corr, coef_t)

    redundant_matrix <- mean_meth_matrix[, correlated_regions]
    non_redundant_matrix <- mean_meth_matrix[, !(colnames(mean_meth_matrix) %in% correlated_regions)]

    print(paste("Number of highly correlated regions", dim(redundant_matrix)[2]))
    print(paste("Number of regions to keep", dim(non_redundant_matrix)[2]))

    if (!is.null(true_file)) {
        true_clone_membership <- read.csv(true_file, sep = "\t", header = TRUE)
        annotation_row <- set_annotation_row(true_clone_membership)
        index_gaps <- set_index_gaps(true_clone_membership)

        print("Finding and removing regions with low variance between clusters")
        low_variance_regions <- find_low_v_regions(non_redundant_matrix, true_clone_membership)

        print(paste("Number of low variance regions", length(low_variance_regions)))
        redundant_matrix <- cbind(redundant_matrix, non_redundant_matrix[, low_variance_regions])
        non_redundant_matrix <- non_redundant_matrix[, !(colnames(non_redundant_matrix) %in% low_variance_regions)]

        redundant_matrix <- redundant_matrix[order(as.integer(true_clone_membership$epigenotype_id)),]
        non_redundant_matrix <- non_redundant_matrix[order(as.integer(true_clone_membership$epigenotype_id)),]
    } else {
        annotation_row <- NULL
        index_gaps <- set_index_gaps(mean_meth_matrix)
    }

    print(paste("Number of redundant regions", dim(redundant_matrix)[2]))
    print(paste("Number of regions to keep", dim(non_redundant_matrix)[2]))

    print("Plotting methylation of redundant regions")
    redundant_plot <- paste0(outdir, "/redundant_regions_plot.png")
    pheatmap(redundant_matrix, cluster_cols = FALSE, cluster_rows = FALSE,
        annotation_row = annotation_row, fontsize = 8,
        main = paste("Mean methylation of redundant regions with threshold", coef_threshold),
        gap_row = index_gaps, fontsize_row = 4, fontsize_col = 6,
        annotation_names_row = FALSE, border_color = NA, show_colnames = TRUE, labels_col = NULL,
        filename = redundant_plot)

    redundant_file <- gzfile(paste0(outdir, "/redundant_regions_meth.tsv.gz"))
    write.csv(redundant_matrix, file = redundant_file, sep = "\t", row.names = TRUE, col.names = TRUE)

    print("Plotting methylation of regions to keep")
    non_redundant_plot <- paste0(outdir, "/regions_to_keep_plot.png")
    pheatmap(non_redundant_matrix, cluster_cols = FALSE, cluster_rows = FALSE,
        annotation_row = annotation_row, fontsize = 8,
        main = paste("Mean methylation of regions to keep with threshold", coef_threshold),
        gap_row = index_gaps, fontsize_row = 4, fontsize_col = 6,
        annotation_names_row = FALSE, border_color = NA, show_colnames = TRUE, labels_col = NULL,
        filename = non_redundant_plot)

    non_redundant_file <- gzfile(paste0(outdir, "/to_keep_meth.tsv.gz"))
    write.csv(non_redundant_matrix, file = non_redundant_file, sep = "\t", row.names = TRUE, col.names = TRUE)

    print("Plotting histogram of Pearson correlation coefficients")
    png(paste0(outdir, "/pearson_corr_coeff_hist.png"))
    hist(as.vector(as.dist(pearson_corr)), breaks = 100)
    dev.off()
}

main(mean_meth_file, true_file, output_directory, coef_t)

