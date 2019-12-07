.load_vis_data <- function(input_CpG_data, input_CN_data_file=NULL, inferred_clusters_file=NULL, true_clusters_file=NULL) {
  # Copy number data
  input_CN_data <- NULL
  if (!is.null(input_CN_data_file)) {
    tmp <- read.csv(input_CN_data_file, sep="\t", header=TRUE, check.names=FALSE)
    input_CN_data <- as.matrix(tmp[,-1])
    rownames(input_CN_data) <- tmp$cell_id
    colnames(input_CN_data) <- sapply(strsplit(colnames(input_CN_data), ":"), function(x) {return(x[1])})
    rm(tmp)
  }

  # inferred clusters from Epiclomal, if available
  inferred_cell_clusters <- NULL
  if (!is.null(inferred_clusters_file)) {
    tmp <- read.csv(inferred_clusters_file, sep="\t", header=TRUE, check.names=FALSE)
    colnames(tmp) <- c("cell_id", "inferred_clusters")
    inferred_cell_clusters <- as.matrix(tmp[,2])
    rownames(inferred_cell_clusters) <- tmp$cell_id
    inferred_cell_clusters <- inferred_cell_clusters[rownames(input_CpG_data),,drop=FALSE]
    inferred_cell_clusters <- as.data.frame(inferred_cell_clusters)
    colnames(inferred_cell_clusters) <- "inferred_clusters"

    if (sum(rownames(input_CpG_data) != rownames(inferred_cell_clusters)) > 0) {
      stop("order of cell IDs 1 doesn't match")
    }
    rm(tmp)
  }

  # True clusters
  true_cell_clusters <- NULL
  if (!is.null(true_clusters_file)) {
    tmp <- read.csv(true_clusters_file, sep="\t", header=TRUE, check.names=FALSE)
    colnames(tmp) <- c("cell_id", "true_clusters")
    true_cell_clusters <- as.matrix(tmp[,2])
    rownames(true_cell_clusters) <- tmp$cell_id
    # MA: reordering inferred_cell_clusters so it matches the order of input_CpG_data
    true_cell_clusters <- true_cell_clusters[rownames(input_CpG_data),,drop=FALSE]
    true_cell_clusters <- as.data.frame(true_cell_clusters)
    colnames(true_cell_clusters) <- "true_clusters"
    if (sum(rownames(input_CpG_data) != rownames(true_cell_clusters)) > 0) {
      stop("order of cell IDs 2 doesn't match")
    }
    rm(tmp)
  }
  return(list("input_CN_data" = input_CN_data, "inferred_cell_clusters" = inferred_cell_clusters, "true_cell_clusters" = true_cell_clusters))
}

.set_index_gaps <- function(data) {
  index <- 1:dim(data)[1]
  index_gaps <- index[!duplicated(data[order(data),])] - 1
  index_gaps <- index_gaps[which(index_gaps != 0)]
  return(index_gaps)
}

.set_annotation_row <- function(true_cell_clusters=NULL, inferred_cell_clusters=NULL, name) {
  # annotating the rows by clusters
  # MA: we want to add the true clusters annotation no matter whether we order by true or predicted
  if (!is.null(true_cell_clusters)) {
    if (!is.null(inferred_cell_clusters)) {
      if (name == "InHouse") {
        true <- true_cell_clusters
        true[true==1] <- "SA501 TNBC"
        true[true==2] <- "SA532 ER+PR-Her2+"
        true[true==3] <- "SA609 TNBC"
        pred <- inferred_cell_clusters
        # this depends of Epiclomal cluster labels
        pred[pred==0] <- "cl. 2"
        pred[pred==7] <- "cl. 1"
        pred[pred==5] <- "cl. 3"
        annotation_row <- cbind(pred, true)
        colnames(annotation_row) <- c("Epiclomal", "Patient")
        annotation_row$Epiclomal <- as.factor(annotation_row$Epiclomal)
        annotation_row$Patient <- as.factor(annotation_row$Patient)
      } else {
        annotation_row <- cbind(inferred_cell_clusters, true_cell_clusters$true_clusters)
        colnames(annotation_row) <- c("inferred clusters", "true clusters")
        annotation_row$`inferred clusters` <- as.factor(annotation_row$`inferred clusters`)
        annotation_row$`true clusters` <- as.factor(annotation_row$`true clusters`)
      }
    } else {
      annotation_row <- true_cell_clusters
      colnames(annotation_row) <- "true clusters"
      annotation_row$`true clusters` <- as.factor(annotation_row$`true clusters`)
    }
  } else {
    annotation_row <- inferred_cell_clusters
    colnames(annotation_row) <- "inferred clusters"
    annotation_row$`inferred clusters` <- as.factor(annotation_row$`inferred clusters`)
  }
  return(annotation_row)
}

.set_plot_data <- function(order_by_true, input_data, true_cell_clusters=NULL, inferred_cell_clusters=NULL) {
  if (order_by_true == 1) {
    data <- input_data[order(as.integer(true_cell_clusters$true_clusters)),]
  } else {
    data <- input_data[order(as.integer(inferred_cell_clusters$inferred_clusters)),]
  }
  return(data)
}

##====================================================================
#' Visualization for Epiclomal
#'
#' @export
#'
#' @param outdir Path to output directory
#' @param input_CpG_data_file Path to methylation data
#' @param input_regions_file Path to region coordinates
#' @param input_CN_data_file Path to copy number data if available
#' @param inferred_clusters_file Path to inferred clusters if available. The third column has posterior probabilities.
#' @param true_clusters_file Path to true clusters if available
#' @param order_by_true If 1, rows are ordered by the true clustering if given, else by the predicted
#' @param name A name for the final plot
#' @param regions_to_plot A file with which regions to plot, for example the flipped_regions.tsv file from the generator.
#' @param use_cache if 1, use cached data when available, if 0, recompute data
#'
#'
#' @importFrom pheatmap pheatmap
#'
#' @examples
#'
#' visualization(outdir, input_CpG_data_file, input_CN_data_file, input_regions_file, inferred_clusters_file, true_clusters_file, order_by_true, name, regions_to_plot, use_cache)
#'
visualization <- function(outdir, input_CpG_data_file, input_regions_file, input_CN_data_file=NULL, inferred_clusters_file=NULL, true_clusters_file=NULL, order_by_true=0, name, regions_to_plot=NULL) {
  dir.create(outdir, showWarnings=FALSE)

  # We can have the following situations:
  # - (NOT YET IMPLEMENTED) neither true_clusters or inferred_clusters are given (e.g. if we want to visualize real data)
  # - only true_clusters are given
  # - only inferred_clusters are given
  # - both true_clusters and inferred_clusters are given, in which case we can order the rows (cells) by true or by inferred

  # Make sure the order_by_true makes sense
  if (!is.null(true_clusters_file) && is.null(inferred_clusters_file)) {
    order_by_true = 1
  }
  if (is.null(true_clusters_file) && !is.null(inferred_clusters_file)) {
    order_by_true = 0
  }

  #======================
  # loading the data
  #======================
  data <- load_data(outdir, input_CpG_data_file, input_regions_file, use_cache=0)

  input_CpG_data <- data$input_CpG_data
  input_regions <- data$input_regions
  mean_meth_matrix <- data$mean_meth_matrix
  rm(data)

  R <- dim(input_regions)[1] ## number of regions
  M <- dim(input_CpG_data)[2] ## number of loci

  data <- .load_vis_data(input_CpG_data, input_CN_data_file=input_CN_data_file, inferred_clusters_file=inferred_clusters_file, true_clusters_file=true_clusters_file)
  input_CN_data <- data$input_CN_data
  inferred_cell_clusters <- data$inferred_cell_clusters
  true_cell_clusters <- data$true_cell_clusters
  rm(data)

  print("number of regions:")
  print(R)
  print("number of loci:")
  print(M)

  #======================
  # pheatmap plots
  #======================
  if (!is.null(true_clusters_file)) {
    if (order_by_true == 1 || is.null(inferred_clusters_file)) {
      index_gaps <- .set_index_gaps(true_cell_clusters)
    }
  }
  if (!is.null(inferred_clusters_file)) {
    if (order_by_true == 0 || is.null(true_clusters_file)) {
      index_gaps <- .set_index_gaps(inferred_cell_clusters)
    }
  }

  annotation_row <- .set_annotation_row(true_cell_clusters=true_cell_clusters, inferred_cell_clusters=inferred_cell_clusters, name=name)

  #======================
  # Plotting methylation data
  #======================

  # annotating the columns by regions for CpG based plots
  reg_id <- unlist(sapply(1:R, function(x){rep(x, (input_regions[x,2]-input_regions[x,1])+1)}))

  annotation_col <- as.matrix(reg_id, nrow=length(colnames(input_CpG_data)))

  rownames(annotation_col) <- colnames(input_CpG_data)
  annotation_col <- as.data.frame(annotation_col)
  colnames(annotation_col) <- "regions"
  annotation_col$regions <- as.factor(annotation_col$regions)

  # adding chromosome labels to the columns

  if (grepl(colnames(input_CpG_data)[1], pattern=":")) {
    show_col_chr_labels <- TRUE
    labels_tmp <- as.vector(sapply(rownames(annotation_col), function(x){unlist(strsplit(x, ":"))[1]}))
    labels_col <- rep(NA, length(labels_tmp))
    labels_col[!duplicated(labels_tmp)] <- labels_tmp[!duplicated(labels_tmp)]
    labels_col[is.na(labels_col)] <- ""
  } else {
    show_col_chr_labels <- FALSE
    labels_col <- NULL
  }

  if (M <= 250) {
    print("Plotting CpG based data")

    data <- .set_plot_data(order_by_true, input_CpG_data, true_cell_clusters, inferred_cell_clusters)

    pheatmap(data, cluster_rows=FALSE, cluster_cols=FALSE, annotation_row=annotation_row, cellwidth=5,
      cellheight=5, fontsize=8,
      main=paste("CpG based methylation data for", name),
      gaps_row=index_gaps, fontsize_row=6, fontsize_col=4,
      annotation_names_row=FALSE, annotation_names_col=FALSE,
      show_colnames=FALSE,
      annotation_col=annotation_col,
      filename=file.path(outdir, paste0(name, "_CpG_based_PLOT.pdf")))

    rm(data)
  }

  if (M > 250) {
    if (R == 1) {
      print("plotting CpG based methylation data for R=1")

      data <- .set_plot_data(order_by_true, input_CpG_data, true_cell_clusters, inferred_cell_clusters)

      pheatmap(data, cluster_cols=FALSE, annotation_row=annotation_row,
        cluster_rows=FALSE,
        fontsize=8, main=paste("CpG based methylation data for", name),
        gaps_row=index_gaps, fontsize_row=8, fontsize_rcol=4,
        annotation_names_row=TRUE,
        show_colnames=show_col_chr_labels,
        labels_col=labels_col,
        filename=file.path(outdir, paste0(name, "_CpG_based_PLOT.pdf")))

      rm(data)
    }

    if (R > 1) {
      if (M < 8000) {
        print("plotting CpG based methylation data for R > 1")

        data <- .set_plot_data(order_by_true, input_CpG_data, true_cell_clusters, inferred_cell_clusters)

        pheatmap(data, cluster_cols=FALSE, annotation_row = annotation_row,
          cluster_rows = FALSE,
          #cellwidth = 5, cellheight = 5,
          fontsize = 8, main = paste0("CpG based methylation data for ", name),
          gaps_row = index_gaps,fontsize_row=8,fontsize_col=6,
          annotation_names_row = FALSE,
          #annotation_colors = ann_colors,
          #annotation_names_col= TRUE,
          show_colnames=show_col_chr_labels,
          labels_col=labels_col,
          #annotation_col=annotation_col,
          annotation_legend = TRUE,
          #legend_breaks = c(0,1),legend_labels = c("unmeth","meth"),
          filename = file.path(outdir, paste0(name,"_CpG_based_PLOT.pdf")))

        rm(data)
      }

      print("plotting region based mean methylation data")

      data <- .set_plot_data(order_by_true, mean_meth_matrix, true_cell_clusters, inferred_cell_clusters)

      if (nrow(data) < 100) {
        fontrow = 8
      } else if (nrow(data) < 200) {
        fontrow = 4
      } else {
        fontrow = 1
      }

      # If we decide to annotate by columns we need to create the object annotation_col for the region based case

      if (show_col_chr_labels) {
        show_col_chr_labels_reg <- TRUE

        CpG_names <- colnames(input_CpG_data)

        labels_tmp <- as.vector(sapply(CpG_names, function(x){unlist(strsplit(x, ":"))[1]}))
        regions_tmp <- labels_tmp[!duplicated(reg_id)]

        labels_col_reg <- rep(NA, length(regions_tmp))
        labels_col_reg[!duplicated(regions_tmp)] <- regions_tmp[!duplicated(regions_tmp)]

        labels_col_reg[is.na(labels_col_reg)] <- ""
      } else {
        show_col_chr_labels_reg <- FALSE
        labels_col_reg <- NULL
      }

      save(data, file=file.path(outdir, paste0(name, "_region_based_data.Rda")))

      if (name == "InHouse") {
        ann_colors = list(Epiclomal = c(`cl. 1`="orange3", `cl. 2`="seagreen", `cl.3`="royalblue"))

        pheatmap(data, cluster_cols=FALSE,
          cluster_rows=FALSE,
          annotation_row=annotation_row,
          annotation_colors=ann_colors,
          fontsize=12,
          main=paste("Mean methylation data for", name),
          gaps_row=index_gaps, fontsize_row=fontrow, fontsize_col=8,
          annotation_names_row=FALSE,
          border_color=NA,
          show_colnames=show_col_chr_labels_reg,
          show_rownames=FALSE,
          labels_col=labels_col_reg,
          filename=file.path(outdir, paste0(name, "_region_based_PLOT.pdf")))
      } else {
        pheatmap(data, cluster_cols=FALSE,
          cluster_rows=FALSE,
          annotation_row=annotation_row,
          fontsize=8,
          main=paste("Region-based mean methylation fraction data for", name),
          gaps_row=index_gaps, fontsize_row=fontrow, fontsize_col=6,
          annotation_names_row=FALSE,
          border_color=NA,
          show_colnames=show_col_chr_labels_reg,
          labels_col=labels_col_reg,
          filename=file.path(outdir, paste0(name, "_region_based_PLOT.pdf")))
      }

      rm(data)

      #=============================
      ### plotting CpG data from a subset of regions
      #=============================

      if (!is.null(regions_to_plot)) {
        regions_for_plot <- scan(file=regions_to_plot, what=integer())
        print("Plotting CpGs for regions")
        print(regions_for_plot)
        input_regions <- input_regions[sort(regions_for_plot),]
        print(input_regions)
        R <- length(regions_for_plot) ## number of regions
        diff_CpG_data <- NULL
        for (r in 1:R) {
          diff_CpG_data <- cbind(diff_CpG_data, input_CpG_data[,c(input_regions[r,1]:input_regions[r,2])])
        }

        M <- dim(diff_CpG_data)[2] ## number of loci

        reg_id <- unlist(sapply(1:R, function(x){rep(x, (input_regions[x,2]-input_regions[x,1])+1)}))
        annotation_col <- as.matrix(reg_id, nrow=length(colnames(diff_CpG_data)))
        rownames(annotation_col) <- colnames(diff_CpG_data)
        annotation_col <- as.data.frame(annotation_col)
        colnames(annotation_col) <- "regions"
        annotation_col$regions <- as.factor(annotation_col$regions)

        data <- .set_plot_data(order_by_true, diff_CpG_data, true_cell_clusters, inferred_cell_clusters)

        pheatmap(data, cluster_rows=FALSE, cluster_cols=FALSE,
          cellwidth=5, cellheight=5,
          fontsize=8, main=paste("CpG-based methylation data for", name),
          gaps_row=index_gaps, fontsize_row=6, fontsize_col=4,
          annotation_names_row=FALSE, annotation_names_col=FALSE,
          show_colnames=FALSE, annotation_col=annotation_col,
          gaps_col=(which(!duplicated(reg_id) == TRUE)[-1]-1),
          filename=file.path(outdir, paste0(name, "_flipped_regions_CpG_based_PLOT.pdf")))

        rm(data)
      }
    }
  }

  #======================
  # Plotting CN data if available
  #======================

  if (!is.null(input_CN_data_file)) {
    ## annotating the columns by chr, but it is currently not working very well
    annotation_col_chr <- as.matrix(colnames(input_CN_data), nrow=length(colnames(test)))
    rownames(annotation_col_chr) <- colnames(input_CN_data)
    annotation_col_chr <- as.data.frame(annotation_col_chr)
    colnames(annotation_col_chr) <- "chr"
    annotation_col_chr$chr <- paste0("chr", annotation_col_chr$chr)
    annotation_col_chr$chr <- as.factor(annotation_col_chr$chr)
    levels(annotation_col_chr$chr) <- paste0("chr", c(1:22, "X"))

    if (dim(input_CN_data)[2] < 250) {
      tmp <- input_CN_data[order(inferred_cell_clusters),]

      pheatmap(tmp, cluster_rows=FALSE, cluster_cols=FALSE, annotation_row=annotation_row,
        cellwidth=6, cellheight=6,
        fontsize=8, main=paste("Copy number data for", name),
        gaps_row=index_gaps, fontsize_row=6, fontsize_col=6, annotation_names_row=FALSE,
        filename=file.path(outdir, paste0(name, "_CN_PLOT.pdf")))

      rm(tmp)
    }

    tmp <- input_CN_data[order(inferred_cell_clusters),]

    pheatmap(tmp,cluster_rows = FALSE,cluster_cols=FALSE, annotation_row = annotation_row,
             fontsize = 8, main = paste0("Copy number data for ", name),
             gaps_row = index_gaps,fontsize_row=6,fontsize_col=4,
             border_color=NA, annotation_names_row = FALSE,show_colnames=FALSE,
             annotation_col=annotation_col_chr,
             #annotation_colors = ann_colors,
             filename = file.path(outdir, paste0(name,"_noLines_CN_PLOT.pdf")))

    rm(tmp)
  }
}
