suppressMessages(library(argparse))
suppressMessages(library(pheatmap))

parser <- ArgumentParser()

parser$add_argument("--mean_methylation_file", default="None", type="character", help="file path contianing mean methylation per region")
parser$add_argument("--true_file", default=NULL, type="character", help="Path to true cluster information")
parser$add_argument("--data_ID", type="character", help="ID of dataset")
parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--coef_threshold", default=0.95, type="double",help="Threshold of Pearson correlation coefficient")
parser$add_argument("--mean_diff_threshold", default=0.05, type="double",help="Threshold of max difference between mean region cluster methylation")


args <- parser$parse_args()

print(args)

data_id <- args$data_ID
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

find_correlated_regions <- function(mean_meth_matrix, coef_t) {
    correlated_regions <- NULL
    for (r1 in 1:(ncol(mean_meth_matrix) - 1)) {
        if (colnames(mean_meth_matrix)[r1] %in% correlated_regions) {
            next
        }
        for (r2 in (r1+1):ncol(mean_meth_matrix)) {
            if (colnames(mean_meth_matrix)[r2] %in% correlated_regions) {
                next
            }

            pearson_corr <- cor(mean_meth_matrix[,r1], mean_meth_matrix[,r2], method = "pearson", use = "na.or.complete")

            if (!is.na(pearson_corr) && ((pearson_corr > coef_t) || (pearson_corr < (-1 * coef_t)))) {
                if (sum(is.na(mean_meth_matrix[r1])) > sum(is.na(mean_meth_matrix[r1]))) {
                    correlated_regions <- append(correlated_regions, colnames(mean_meth_matrix)[r1])
                    next
                } else {
                    correlated_regions <- append(correlated_regions, colnames(mean_meth_matrix)[r2])
                }
            }
        }
    }
    return(correlated_regions)
}

find_cluster_means <- function(data, true_clusters) {
    data <- t(data)
    epi_means <- NULL
    for (epi in unique(true_clusters$epigenotype_id)) {
        epi_means <- cbind(epi_means, rowMeans(data[,true_clusters$epigenotype_id == epi], na.rm = TRUE))
    }
    epi_means <- as.data.frame(epi_means)
    epi_means$max_diff <- apply(epi_means, 1, function(x) diff(range(x, na.rm = TRUE)))

    return(epi_means)
}

main <- function(mean_meth_file, true_file, output_directory, coef_threshold) {
    mean_meth_matrix <- load_mean_meth(mean_meth_file)
    print(paste("Number of cells", dim(mean_meth_matrix)[1]))
    print(paste("Number of regions", dim(mean_meth_matrix)[2]))

    print("Finding and removing highly correlated regions")
    correlated_regions <- find_correlated_regions(mean_meth_matrix, coef_t)

    redundant_matrix <- mean_meth_matrix[, colnames(mean_meth_matrix) %in% correlated_regions]
    non_redundant_matrix <- mean_meth_matrix[, !(colnames(mean_meth_matrix) %in% correlated_regions)]

    print(paste("Number of highly correlated regions", dim(redundant_matrix)[2]))
    print(paste("Number of regions to keep", dim(non_redundant_matrix)[2]))

    if (!is.null(true_file)) {
        true_clone_membership <- read.csv(true_file, sep = "\t", header = TRUE)
        annotation_row <- set_annotation_row(true_clone_membership)
        index_gaps <- set_index_gaps(true_clone_membership)

        redundant_matrix <- redundant_matrix[order(as.integer(true_clone_membership$epigenotype_id)),]
        non_redundant_matrix <- non_redundant_matrix[order(as.integer(true_clone_membership$epigenotype_id)),]
    } else {
        annotation_row <- NULL
        index_gaps <- set_index_gaps(mean_meth_matrix)
    }

    redundant_file <- gzfile(paste0(outdir, "/redundant_regions_meth_", data_id, ".tsv.gz"))
    write.csv(redundant_matrix, file = redundant_file, sep = "\t", row.names = TRUE, col.names = TRUE)

    non_redundant_file <- gzfile(paste0(outdir, "/to_keep_meth_", data_id, ".tsv.gz"))
    write.csv(non_redundant_matrix, file = non_redundant_file, sep = "\t", row.names = TRUE, col.names = TRUE)

    regions_file <- paste0(outdir, "/filtered_regions_", data_id, ".tsv")
    write.table(colnames(non_redundant_matrix), file = regions_file, sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

    print("Plotting methylation of redundant regions")
    redundant_plot <- paste0(outdir, "/redundant_regions_plot_", data_id, ".png")
    pheatmap(redundant_matrix, cluster_cols = FALSE, cluster_rows = FALSE,
        annotation_row = annotation_row, fontsize = 8,
        main = "Mean methylation of redundant regions",
        gap_row = index_gaps, fontsize_row = 4, fontsize_col = 6,
        annotation_names_row = FALSE, border_color = NA, show_colnames = TRUE, labels_col = NULL,
        filename = redundant_plot)

    print("Plotting methylation of regions to keep")
    non_redundant_plot <- paste0(outdir, "/regions_to_keep_plot_", data_id, ".png")
    pheatmap(non_redundant_matrix, cluster_cols = FALSE, cluster_rows = FALSE,
        annotation_row = annotation_row, fontsize = 8,
        main = "Mean methylation of regions to keep",
        gap_row = index_gaps, fontsize_row = 4, fontsize_col = 6,
        annotation_names_row = FALSE, border_color = NA, show_colnames = TRUE, labels_col = NULL,
        filename = non_redundant_plot)

    if (!is.null(true_file)) {

        print("Finding and removing regions with low variance between clusters")
        cluster_means <- find_cluster_means(mean_meth_matrix, true_clone_membership)
        low_variance_regions <- rownames(cluster_means[cluster_means$max_diff < 0.05,])

        print(paste("Number of low variance regions", length(low_variance_regions)))
        redundant_matrix <- mean_meth_matrix[, colnames(mean_meth_matrix) %in% low_variance_regions]
        non_redundant_matrix <- mean_meth_matrix[, !(colnames(mean_meth_matrix) %in% low_variance_regions)]

        print("Plotting histogram of max difference in methylation between clusters by region")
        png(paste0(outdir, "/mean_diff_hist_", data_id, ".png"))
        hist(cluster_means$max_diff, breaks = 100)
        dev.off()

        redundant_matrix <- redundant_matrix[order(as.integer(true_clone_membership$epigenotype_id)),]
        non_redundant_matrix <- non_redundant_matrix[order(as.integer(true_clone_membership$epigenotype_id)),]

        if (sum(!is.na(redundant_matrix)) > 0) {
            print("Plotting methylation of low variance regions")
            redundant_plot <- paste0(outdir, "/low_variance_regions_plot_", data_id, ".png")
            pheatmap(redundant_matrix, cluster_cols = FALSE, cluster_rows = FALSE,
                annotation_row = annotation_row, fontsize = 8,
                main = "Mean methylation of low variance regions",
                gap_row = index_gaps, fontsize_row = 4, fontsize_col = 6,
                annotation_names_row = FALSE, border_color = NA, show_colnames = TRUE, labels_col = NULL,
                filename = redundant_plot)
        }

        if (sum(!is.na(non_redundant_matrix)) > 0) {
            print("Plotting methylation of variant regions")
            non_redundant_plot <- paste0(outdir, "/variant_regions_plot_", data_id, ".png")
            pheatmap(non_redundant_matrix, cluster_cols = FALSE, cluster_rows = FALSE,
                annotation_row = annotation_row, fontsize = 8,
                main = "Mean methylation of variant regions",
                gap_row = index_gaps, fontsize_row = 4, fontsize_col = 6,
                annotation_names_row = FALSE, border_color = NA, show_colnames = TRUE, labels_col = NULL,
                filename = non_redundant_plot)
        }
    }
    # print("Plotting histogram of Pearson correlation coefficients")
    # png(paste0(outdir, "/pearson_corr_coeff_hist_", data_id, ".png"))
    # hist(as.vector(as.dist(pearson_corr)), breaks = 100)
    # dev.off()
}

main(mean_meth_file, true_file, output_directory, coef_t)

