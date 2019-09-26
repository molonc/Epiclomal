##====================================================================
# REpiclomal v1.0
#
# Camila P.E. de Souza
# Department of Statistical and Actuarial Sciences, University of Western Ontario
#
# Mirela Andronescu
# Department of Molecular Oncology, BC Cancer Reseach Centre
#
#
##====================================================================

.extract_mean_meth_per_cell <- function(cell_data,region_coord) {
  mean_meth <- apply(region_coord,1,function(x) {mean(cell_data[x[1]:x[2]],na.rm=TRUE)})
  mean_meth[is.na(mean_meth)] <- NA
  return(mean_meth)
}

##====================================================================
#' Load methylation data
#'
#' @export
#'
#' @param input_CpG_data_file A file containing methylation data
#' @param input_regions_file A file containing regions for region based clustering
#' @param use_cache if 1, use cached data when available, if 0, recompute data
#'
#' @return A list containing the parsed methylation data as a matrix (input_CpG_data),
#' the input regions (input_regions), and the mean methylation per region as a matrix (mean_meth_matrix)
#'
#' @examples
#' data <- load_data(input_CpG_data_file, input_regions_file)
#'
#' input_CpG_data <- data$input_CpG_data
#' input_regions <- data$input_regions
#' mean_meth_matrix <- data$mean_meth_matrix
#' rm(data)
#'
load_data <- function(input_CpG_data_file, input_regions_file, use_cache) {
  cached_data <- gsub(".tsv.gz", ".RDa.gz", input_CpG_data_file)
  if (file.exists(cached_data) & use_cache) {
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

    mean_meth_matrix <- t(apply(input_CpG_data,1, .extract_mean_meth_per_cell, region_coord=input_regions))

    save(input_CpG_data, input_regions, mean_meth_matrix, file = cached_data, compress = "gzip")
  }
  return(list("input_CpG_data" = input_CpG_data, "input_regions" = input_regions, "mean_meth_matrix" = mean_meth_matrix))
}

##====================================================================#' DensityCut caller
#'
#' @export
#'
#' @param input_CpG_data A matrix of methylation data
#' @param mean_meth_matrix A matrix of mean methylation per region
#' @param R An integer, number of regions
#' @param max_PC An integer, max number of PC components (default=20)
#' @param maxit An integer, max number of iterations for DensityCut (default=100)
#' @param impute if 1, fill in NAs with imputed data, if 0, do nothing
#' @param use_cache if 1, use cached data when available, if 0, recompute data
#' @param outdir file path to output directory
#'
#'
#' @return Predicted clusters from DensityCut algorithm
#'
#' @importFrom pcaMethods pca
#' @importFrom pcaMethods scores
#' @importFrom densitycut DensityCut
#' @importFrom pheatmap pheatmap
#'
#' @examples
#' data <- load_data(input_CpG_data_file, input_regions_file)
#'
#' input_CpG_data <- data$input_CpG_data
#' input_regions <- data$input_regions
#'
#' R <- dim(input_regions)[1] ## number of regions
#' maxpc <- min(20, R)
#' outdir <- args$outdir
#'
#' densitycut.clust(input_CpG_data, mean_meth_matrix, R, max_PC, maxit=100, impute=0, use_cache=1, outdir)
#'
densitycut.clust <- function(input_CpG_data, mean_meth_matrix, R, max_PC=20, maxit=100, impute, use_cache, outdir) {
  if (R == 1) {
    stop("One region, cannot do Region-based densityCut")
  }
  if (R > 1) {
    print("More than one region, region based densityCut")
    if (impute == 1) {
      imputed_file <- paste0(outdir, "/region_based_imputed.RDa.gz")
      if (file.exists(imputed_file) & use_cache) {
        print ("Reading the imputed file")
        load(imputed_file)
        print (" ... done.")
      } else {
        # to remove
        #input <- input[1:100,1:5]

        # replace with average values, for each col
        print("Per region, replacing NAs with average values")
        mean_meth_matrix <- .impute_means(mean_meth_matrix)
        print(" ... done.")

        # eliminate the empty rows (features)
        mean_meth_matrix <- mean_meth_matrix[ rowSums(mean_meth_matrix)!=0, ]
        save(mean_meth_matrix, file = imputed_file, compress = "gzip")
      }
    }

    max_comp <- min(max_PC, R)

    print("number of PC components:")
    print(max_comp)

    t <- try(pca(mean_meth_matrix, method="nipals", nPcs=max_comp))
    if ("try-error" %in% class(t)) {
      stop("Stop! At least one cell has NO data across all regions, can't do PCA")
    } else {
      pc <- t
    }
    pc_scores <- scores(pc)

    cluster.out <- DensityCut(pc_scores, maxit=maxit)
    checking_warning <- capture.output(cluster.out)

    if (checking_warning[1] == "WARNING! not converged ") {
      stop("densitycut didn't converge, not saving results")
    } else {
      print("densitycut converged!")

      possible_clusters <- cbind(rownames(input_CpG_data), cluster.out$cluster)
      possible_clusters <- as.data.frame(possible_clusters)
      colnames(possible_clusters) <- c("cell_id", "DensityCut")

      ofile <- paste0(outdir, "/DensityCut_clusters_Region_based_maxPC_", max_comp, ".tsv")
      write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)
      system(paste("gzip --force", ofile))

      print("Plotting DensityCut clusters")

      inf_clustrs_order <- order(as.integer(cluster.out$cluster))

      annotation_row <- as.data.frame(cluster.out$cluster)
      colnames(annotation_row) <- "inferred clusters"
      annotation_row$`inferred clusters` <- as.factor(annotation_row$`inferred clusters`)
      rownames(annotation_row) <- rownames(mean_meth_matrix)

      plot_file <- paste0(outdir, "/DensityCut_PLOT.pdf")
      if (!file.exists(plot_file)) {
        pheatmap(mean_meth_matrix[inf_clustrs_order,], clusters_row=FALSE, cluster_cols=FALSE,
          cellwidth=8, cellheight=8, fontsize=8,
          annotation_row=annotation_row,
          main = "DensityCut",
          show_colnames=FALSE,
          annotation_names_row=FALSE,
          filename=plot_file
        )
      }
      return(possible_clusters)
    }
  }
}

##====================================================================
#' Euclidean Clustering
#'
#' @export
#'
#' @param input_CpG_data A matrix of methylation data
#' @param mean_meth_matrix A matrix of mean methylation per region
#' @param R An integer, number of regions
#' @param Max_K An integer, maximum number of clusters to be considered when cutting the tree
#' @param index_type "ch" or "gap", Index to be used to choose number of clusters
#' @param impute if 1, impute with the average per region/locus, if it is 0 do nothing
#' @param use_cache if 1, use cached data when available, if 0, recompute data
#' @param outdir file path to output directory
#'
#'
#' @return Predicted clusters from Euclidean clustering algorithm
#'
#' @importFrom NbClust NbClust
#' @importFrom pheatmap pheatmap
#'
#' @examples
#' data <- load_data(input_CpG_data_file, input_regions_file)
#'
#' input_CpG_data <- data$input_CpG_data
#' input_regions <- data$input_regions
#'
#' R <- dim(input_regions)[1] ## number of regions
#' Max_K <- args$Max_K
#' outdir <- args$outdir
#'
#' euclidean.clust(input_CpG_data, mean_meth_matrix, R, Max_K, index_type="ch", impute=0, use_cache=1, outdir)
#'
euclidean.clust <- function(input_CpG_data, mean_meth_matrix, R, Max_K, index_type="ch", impute, use_cache, outdir) {
  if (R == 1) {
    if (impute == 1) {
      stop("Imputing option not implemented for this case")
    }

    print("One region, CpG based hiearchical clustering")

    pairwisedist_file <- paste0(outdir, "/pairwisedist.RDa.gz")
    if (file.exists(pairwisedist_file) & use_cache) {
      print("Loading previously calculated pairwise dist")
      load(pairwisedist_file)
      print("... done.")
    } else {
      pairwisedist <- dist(input_CpG_data, method="euclidean")
      save(pairwisedist, file=pairwisedist_file, compress="gzip")
    }

    if (sum(is.na(pairwisedist)) == 0) {
      hclust_CpG_crash <- 0
      hcluster <- hclust(pairwisedist, method="complete")
      mycl <- cutree(hcluster, k=1:Max_K)
      possible_clusters <- cbind(rownames(input_CpG_data), mycl)
      possible_clusters <- as.data.frame(possible_clusters)
      colnames(possible_clusters) <- c("cell_id", paste0("EuclideanClust_cpg_num_clusters_", 1:Max_K))

      t <- try(NbClust(input_CpG_data, diss=pairwisedist, distance=NULL, min.nc=2, max.nc=Max_K, method="complete", index="cindex"))
      if ("try-error" %in% class(t)) {
        print("Can't find a best partition")
        error_ch_index <- 1
        hcluster_Nb <- t
      }

      write.table(error_ch_index, file=paste0(outdir, "/EuclideanClust_bestpartition_crash.tsv"), row.names=FALSE, col.names=FALSE)

      if (error_ch_index == 0) {
        best_cluster <- hcluster_Nb$Best.partition
        possible_clusters <- cbind(possible_clusters, best_cluster)
        colnames(possible_clusters) <- c(colnames(possible_clusters)[1:(dim(possible_clusters)[2]-1)], paste0("EuclideanClust_best_cluster_", hcluster_Nb$Best.nc[1]))
      }

      ofile <- paste0(outdir, "/EuclideanClust_clusters_CpG_based_maxk_", Max_K, ".tsv")
      write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)
      system(paste("gzip --force", ofile))

      ofile <- paste0(outdir, "/EuclideanClust_cell_order_CpG_based_maxk_", Max_K, ".tsv")
      write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
      system(paste("gzip --force", ofile))
    } else {
      print("Some pairs of cells have no CpG with data in common")
      hclust_CpG_crash <- 1
    }
    write.table(hclust_CpG_crash, file=paste0(outdir, "/EuclideanClust_crash.tsv"), row.names=FALSE, col.names=FALSE)
    return(possible_clusters)
  }

  if (R > 1) {
    print("More than one region, region based hiearchical clustering")
    if (impute == 1) {
      imputed_file <- paste0(outdir, "/EuclideanClust_mean_meth_matrix.Rda.gz")
      if (file.exists(imputed_file) & use_cache) {
        print("Reading the imputed cache file")
        load(imputed_file)
        print("... done.")
      } else {
        print("Per region, replacing NAs with average values")
        mean_meth_matrix <- .impute_means(mean_meth_matrix)
        print("... done.")
        mean_meth_matrix <- mean_meth_matrix[rowSums(mean_meth_matrix) != 0, ]
        save(mean_meth_matrix, file=imputed_file, compress="gzip")
      }
    }
    pairwisedist_region_file <- paste0(outdir, "/pairwisedist_region.RDa.gz")
    if (file.exists(pairwisedist_region_file) & use_cache) {
      print("loading pairwise Euclidean distances")
      load(pairwisedist_region_file)
      print("... done.")
    } else {
      print("Doing hclust in (dis)similarity matrix")
      dist_region <- dist(mean_meth_matrix, method="euclidean")
      pairwisedist_region <- dist(dist_region, method="euclidean")
      save(dist_region, pairwisedist_region, file=pairwisedist_region_file, compress="gzip")
      print("... done.")
    }

    if (sum(is.na(pairwisedist_region)) == 0) {
      hclust_region_crash <- 0
      hcluster <- hclust(pairwisedist_region, method="complete")

      print("Plotting Euclidean clusters")
      plot_file <- paste0(outdir, "/Region_based_EuclideanClust_PLOT.pdf")
      if(!file.exists(plot_file)) {
        pheatmap(as.matrix(dist_region), cluster_rows=TRUE, cluster_cols=TRUE,
          cellwidth=8, cellheight=8, fontsize=8,
          clustering_distance_rows="euclidean",
          clustering_method="complete",
          main="Region-based EuclideanClust",
          filename=plot_file
        )
      }

      mycl <- cutree(hcluster, k=1:Max_K)

      possible_clusters <- cbind(rownames(input_CpG_data), mycl)
      possible_clusters <- as.data.frame(possible_clusters)
      colnames(possible_clusters) <- c("cell_id", paste0("EuclideanClust_region_num_clusters_", 1:Max_K))

      t <- try(NbClust(as.matrix(dist_region), diss=pairwisedist_region, distance=NULL, min.nc=1, max.nc=Max_K, method="complete", index=index_type))
      if ("try-error" %in% class(t)) {
        print("can't use ch or gap index")
        error_ch_index <- 1
      } else {
        error_ch_index <- 0
        hcluster_Nb <- t
      }

      write.table(error_ch_index, file=paste0(outdir, "/EuclideanClust_bestpartition_crash.tsv"), row.names=FALSE, col.names=FALSE)

      if (error_ch_index == 0) {
        best_cluster <- hcluster_Nb$Best.partition

        possible_clusters <- cbind(possible_clusters, best_cluster)
        colnames(possible_clusters) <- c(colnames(possible_clusters)[1:(dim(possible_clusters)[2]-1)], paste0("EuclideanClust_best_cluster_", hcluster_Nb$Best.nc[1]))
      }

      ofile <- paste0(outdir, "/EuclideanClust_clusters_region_based_maxk_", Max_K, ".tsv")
      write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)
      system(paste("gzip --force", ofile))

      ofile <- paste0(outdir, "/EuclideanClust_cell_order_region_based_maxk_", Max_K, ".tsv")
      write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
      system(paste("gzip --force", ofile))
    } else {
      print("some pairs of cells have no region with data in common")
      hclust_region_crash <- 1
    }
    write.table(hclust_region_crash, file=paste0(outdir, "/EuclideanClust_crash.tsv"), row.names=FALSE, col.names=FALSE)
  }
  return(possible_clusters)
}

##====================================================================
#' Hamming Clustering
#'
#' @export
#'
#' @param input_CpG_data A matrix of methylation data
#' @param Max_K An integer, maximum number of clusters to be considered when cutting the tree
#' @param index_type "ch" or "gap", Index to be used to choose number of clusters. Default "ch"
#' @param impute if 1, impute with the average per region/locus, if it is 0 do nothing
#' @param use_cache if 1, use cached data when available, if 0, recompute data
#' @param outdir file path to output directory
#'
#'
#' @return Predicted clusters from Euclidean clustering algorithm
#'
#' @importFrom NbClust NbClust
#' @importFrom pheatmap pheatmap
#'
#' @examples
#' data <- load_data(input_CpG_data_file, input_regions_file)
#'
#' input_CpG_data <- data$input_CpG_data

#' Max_K <- args$Max_K
#' outdir <- args$outdir
#'
#' hamming.clust(input_CpG_data, Max_K, index_type="ch", impute=0, use_cache=1, outdir)
#'
hamming.clust <- function(input_CpG_data, Max_K, index_type="ch", impute, use_cache, outdir) {
  print("Tony's approach - CpG based clustering (HammingClust)")

  if (impute == 1) {
    imputed_file <- paste0(outdir, "/CpG_based_imputed.RDa.gz")
    if (file.exists(imputed_file) & use_cache) {
      print("Reading the imputed file")
      load(imputed_file)
      print("... done.")
    } else {
      # replace with average values for each col
      print("Per locus, replacing NAs with median values")
      input_CpG_data <- .impute_medians(input_CpG_data)
      print("... done.")

      # eliminate the empty rows (features)
      input_CpG_data <- input_CpG_data[rowSums(input_CpG_data) != 0, ]
      save(input_CpG_data, file=imputed_file, compress = "gzip")
    }
  }
  dist_PBAL_file <- paste0(outdir, "/dist_PBAL.RDa.gz")
  if (file.exists(dist_PBAL_file) & use_cache) {
    print("Loading PBAL distance matrix from file")
    load(dist_PBAL_file)
    print("... done.")
  } else {
    print("Computing PBAL distance matrix")
    dist_PBAL <- as.dist(.dist_PBAL(d = input_CpG_data))
    print("Computing pairwise PBAL distance matrix")
    diss_matrix <- dist(dist_PBAL, method="euclidean")
    print("Done computing, saving to file")
    save(dist_PBAL, diss_matrix, file=dist_PBAL_file, compress="gzip")
  }

  if (sum(is.na(diss_matrix)) == 0) {
    PBAL_crash <- 0

    print("Plotting Hamming clusters")
    plot_file <- paste0(outdir, "/HammingClust_PLOT.pdf")
    if (!file.exists(plot_file)) {
      pheatmap(dist_PBAL, cluster_rows=TRUE, cluster_cols=TRUE,
        cellwidth=8, cellheight=8, fontsize=8,
        clustering_distance_rows="euclidean",
        cllustering_method="ward.D2",
        main="HammingClust",
        filename=plot_file
      )
    }

    hcluster <- hclust(diss_matrix, method="ward.D2")
    mycl <- cutree(hcluster, k=1:Max_K)

    possible_clusters <- cbind(rownames(input_CpG_data), mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id", paste0("HammingClust_num_clusters_", 1:Max_K))

    t <- try(NbClust(dist_PBAL, diss=diss_matrix, distance=NULL, min.nc=1, max.nc=Max_K, method="ward.D2", index=index_type))
    if ("try-error" %in% class(t)) {
      print("Can't find best partition")
      error_ch_index <- 1
    } else {
      error_ch_index <- 0
      hcluster_Nb <- t
    }

    write.table(error_ch_index, file=paste0(outdir, "/HammingClust_bestpartition_crash.tsv"), row.names=FALSE, col.names=FALSE)

    if (error_ch_index == 0) {
      best_cluster <- hcluster_Nb$Best.partition

      possible_clusters <- cbind(possible_clusters, best_cluster)
      colnames(possible_clusters) <- c(colnames(possible_clusters)[1:(dim(possible_clusters)[2]-1)], paste0("HammingClust_best_cluster_", hcluster_Nb$Best.nc[1]))
    }

    ofile <- paste0(outdir, "/HammingClust_clusters_CpG_based_maxk_", Max_K, ".tsv")
    write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)
    system(paste("gzip --force", ofile))

    ofile <- paste0(outdir, "/HammingClust_cell_order_CpG_based_maxk_", Max_K, ".tsv")
    write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste("gzip --force", ofile))
  } else {
    PBAL_crash <- 1
  }

  write.table(PBAL_crash, file=paste0(outdir, "/HammingClust_crash.tsv"), row.names=FALSE, col.names=FALSE)
  return(possible_clusters)
}

##====================================================================
#' Pearson Clustering
#'
#' @export
#'
#' @param input_CpG_data A matrix of methylation data
#' @param Max_K An integer, maximum number of clusters to be considered when cutting the tree
#' @param index_type "ch" or "gap", Index to be used to choose number of clusters. Default "ch"
#' @param impute if 1, impute with the average per region/locus, if it is 0 do nothing
#' @param use_cache if 1, use cached data when available, if 0, recompute data
#' @param outdir file path to output directory
#'
#'
#' @return Predicted clusters from Euclidean clustering algorithm
#'
#' @importFrom NbClust NbClust
#' @importFrom pheatmap pheatmap
#'
#' @examples
#' data <- load_data(input_CpG_data_file, input_regions_file)
#'
#' input_CpG_data <- data$input_CpG_data

#' Max_K <- args$Max_K
#' outdir <- args$outdir
#'
#' pearson.clust(input_CpG_data, Max_K, index_type="ch", impute=0, use_cache=1, outdir)
#'
pearson.clust <- function(input_CpG_data, Max_K, index_type="ch", impute, use_cache, outdir) {
  if (impute == 1) {
    imputed_file <- paste0(outdir, "/CpG_based_imputed.RDa.gz")
    if (file.exists(imputed_file) & use_cache) {
      print("Loading the imputed file")
      load(imputed_file)
      print("... done.")
    } else {
      print("Per region, replacing NAs with median values")
      input_CpG_data <- .impute_medians(input_CpG_data)
      print("... done.")

      input_CpG_data <- input_CpG_data[rowSums(input_CpG_data) != 0, ]
      save(input_CpG_data, file=imputed_file, compress="gzip")
    }
  }

  dist_Pearson_file <- paste0(outdir, "/dist_Pearson.RDa.gz")
  if (file.exists(dist_Pearson_file) & use_cache) {
    print("Loading Pearson distances from file")
    load(dist_Pearson_file)
    print("... done.")
  } else {
    print("computing Pearson distances")
    dist_Pearson <- cor(x=t(input_CpG_data), method="pearson", use="pairwise.complete.obs")
    print("scTrio's apprioach - CpG based clustering")
    diss_matrix <- 1 - cor(x=dist_Pearson, method="pearson", use="pairwise.complete.obs")
    print("done computing, saving to file")
    save(dist_Pearson, diss_matrix, file=dist_Pearson_file, compress="gzip")
  }

  if (sum(is.na(diss_matrix)) == 0) {
    Pearson_crash <- 0

    hcluster <- hclust(as.dist(diss_matrix), method = "ward.D2")

    print("Plotting Pearson Clusters")
    plot_file <- paste0(outdir, "/PearsonClust_PLOT.pdf")
    if (!file.exists(plot_file)) {
      pheatmap(dist_Pearson, cluster_rows=TRUE, cluster_cols=TRUE,
        cellwidth=8, cellheight=8, fontsize=8,
        clustering_distance_cols=as.dist(diss_matrix),
        clustering_distance_rows=as.dist(diss_matrix),
        clustering_method="ward.D2",
        main="Pearson corr. approach",
        filename=plot_file
      )
    }

    mycl <- cutree(hcluster, k=1:Max_K)

    possible_clusters <- cbind(rownames(input_CpG_data), mycl)
    possible_clusters <- as.data.frame(possible_clusters)
    colnames(possible_clusters) <- c("cell_id", paste0("PearsonClust_num_clusters_", 1:Max_K))

    t <- try(NbClust(dist_Pearson, diss=as.dist(diss_matrix), distance=NULL, min.nc=1, max.nc=Max_K, method="ward.D2", index=index_type))
    if ("try-error" %in% class(t)) {
      print("can't find best partition")
      error_ch_index <- 1
    } else {
      error_ch_index <- 0
      hcluster_Nb <- t
    }

    write.table(error_ch_index, file=paste0(outdir, "/PearsonClust_bestpartition_crash.tsv"), row.names=FALSE, col.names=FALSE)

    if (error_ch_index == 0) {
      best_cluster <- hcluster_Nb$Best.partition

      possible_clusters <- cbind(possible_clusters, best_cluster)
      colnames(possible_clusters) <- c(colnames(possible_clusters)[1:(dim(possible_clusters)[2]-1)], paste0("PearsonClust_best_cluster_", hcluster_Nb$Best.nc[1]))
    }

    ofile <- paste0(outdir, "/PearsonClust_clusters_CpG_based_maxk_", Max_K, ".tsv")
    write.table(possible_clusters, file=ofile, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)
    system(paste("gzip --force", ofile))

    ofile <- paste0(outdir, "/PearsonClust_cell_order_CpG_based_maxk_", Max_K, ".tsv")
    write.table(hcluster$order, file=ofile, sep="\t", col.names=FALSE, quote=FALSE)
    system(paste("gzip --force", ofile))
  } else {
    Pearson_crash <- 1
  }
  write.table(Pearson_crash, file=paste0(outdir, "/PearsonClust_crash.tsv"), row.names=FALSE, col.names=FALSE)
  return(possible_clusters)
}



