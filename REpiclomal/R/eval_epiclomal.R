######################################

### auxiliary functions
### https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve

elbow_finder <- function(x_values, y_values) {
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)

  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }

  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]

  print("In elbow_finder, x, y, distances")
  print(x_values)
  print(y_values)
  print(distances)
  print(paste0("Max dist ", x_max_dist))

  #return(distances)
  #return(c(x_max_dist, y_max_dist))
  return(c(x_max_dist))
}

extract_yaml <- function(input_data, data_field) {
  values <- sapply(input_data, '[[', data_field)
  values[sapply(values, is.null)] <- NA
  values <- as.numeric(sapply(values, unlist))
  return(values)
}

##====================================================================
#' Evaluation function for epiclomal results
#'
#' @export
#'
#' @param input input directory of epiclomal results
#' @param outdir output directory of evaluation results
#' @param model A name for the model
#' @param flag flag
#' @param criterion criterion
#' @param GAIN_THRESHOLD gain threshold
#'
#'
#' @return list of best row and table best. as best_row and table_best
#'
#' @importFrom yaml yaml.load_file
#'
#' @examples
#' run_eval (input, "all", "DIC_LINE_ELBOW", "0.05_-100")

evaluate.epiclomal <- function (input, outdir, model, flag, criterion, GAIN_THRESHOLD) {
  # criterion can be "lower_bound" or "log_posterior"
  outdir <- file.path(outdir, paste0(criterion, "_gainthr", GAIN_THRESHOLD))
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

  # input can actually be a list of directories. Then look through all of them and compute the measure
  directories <- Sys.glob(input)
  files <- unlist(lapply(directories, function(x) return(list.files(x, recursive=TRUE, pattern="params.yaml", full.names=TRUE))))
  print(paste("number of files:", length(files)))
  run <- unlist(lapply(directories, function(x) return(list.files(x, recursive=TRUE, pattern="cluster_posteriors.tsv.gz", full.names=TRUE))))
  print(paste("length of runs:", length(run)))
  lines <- lapply(files, yaml.load_file)
  print(paste("length of lines:", length(lines)))

  converged <- as.integer(as.logical(sapply(lines, '[[', "converged")))
  if (criterion == "DIC_measure" || criterion == "DIC_LINE_ELBOW") {
    measure <- "DIC_measure"
  } else {
    measure <- criterion
  }

  score <- extract_yaml(lines, measure)
  cpu_time <- extract_yaml(lines, "CPU_time_seconds")
  memory <- extract_yaml(lines, "Max_memory_MB")
  all_vmeasure <- extract_yaml(lines, "Vmeasure")
  nclusters_pred <- extract_yaml(lines, "nclusters")
  # for clone_prev_MAE, I may also have slsbulk_clone_prev_MAE
  clone_prev_MAE <- extract_yaml(lines, "clone_prev_MAE")
  clone_prev_MSE <- extract_yaml(lines, "clone_prev_MSE")
  slsbulk_vmeasure <- extract_yaml(lines, "slsbulk_vmeasure")
  slsbulk_clone_prev_MAE <- extract_yaml(lines, "slsbulk_clone_prev_MAE")
  slsbulk_clone_prev_MSE <- extract_yaml(lines, "slsbulk_clone_prev_MSE")
  uncertainty <- extract_yaml(lines, "uncertainty_true_positive_rate")

  table_all <- data.frame(converged, score, run, cpu_time, memory, nclusters_pred, all_vmeasure, clone_prev_MAE, clone_prev_MSE, slsbulk_clone_prev_MSE, slsbulk_vmeasure, slsbulk_clone_prev_MAE, uncertainty)

  dfile <- file.path(outdir, paste0(flag, "_results_allruns_", model, ".tsv"))
  write.table(table_all, file = dfile, quote = FALSE, sep = "\t", row.names = FALSE)

  ntotal <- nrow(table_all)
  nconv <- sum(table_all$converged)

  if (criterion == "DIC_measure" || criterion == "DIC_LINE_ELBOW") {
    table_per_cluster = data.frame()
    for (k in sort(unique(table_all$nclusters_pred))) {
      rows <- table_all[which(table_all$nclusters_pred==k),]
      minofk_row <- rows[which.min(rows$score),]
      table_per_cluster = rbind(table_per_cluster, minofk_row)
    }
    # Now measure the relative gain
    print("Table per cluster")
    print(table_per_cluster)

    dfile <- file.path(outdir, paste0(flag, "_results_perclusterruns_", model, ".tsv"))
      write.table(table_per_cluster, file=dfile, quote=FALSE, sep="\t", row.names=FALSE)

      if (sum(!is.na(table_per_cluster$all_vmeasure))) {
      # also plot the v-measure versus the number of predicted clusters to see if we get better V-measure when we choose a different number of clusters
      pdf(file.path(outdir, "Vmeasure_vs_nclusters.pdf"), height=7, width=9)
      x <- table_per_cluster$nclusters_pred
      y <- table_per_cluster$all_vmeasure

      # type='o' means it plots both points and lines overplotted
      matplot(x, y, lty=1, type="o", lwd=c(4), col=c(4), pch=19,
        ylab="V-measure for run with best score",
        xlab="Number of clusters",
        cex.axis=1.2,cex.lab=1.2)
      dev.off()
    }

    print(paste("Number of rows in table per cluster is", nrow(table_per_cluster)))

    if (nrow(table_per_cluster) == 1) {
      best_score <- table_per_cluster[1,2]
      best_row <- table_per_cluster[1,]
    } else {
      if (criterion == "DIC_measure") {
        for (k in 1:(nrow(table_per_cluster)-1)) {
          gain <- (table_per_cluster[k,2] - table_per_cluster[k+1,2])/table_per_cluster[k,2]
          print(paste("Gain", gain))
          if (gain < GAIN_THRESHOLD) {
            best_score <- table_per_cluster[k,2]
            best_row <- table_per_cluster[k,]
          } else {
            best_score <- table_per_cluster[k+1, 2]
            best_row <- table_per_cluster[k+1,]
          }
        }
      }
      if (criterion == "DIC_LINE_ELBOW") {
        print("Doing DIC LINE ELBOW")
        THRESHOLD <- as.numeric(unlist(strsplit(GAIN_THRESHOLD, split="_")))
        x_elbow <- c(table_per_cluster[1,c("nclusters_pred")])
        y_elbow <- c(table_per_cluster[1,c("score")])
        gain_vector <- c(0)
        for (k in 1:(nrow(table_per_cluster)-1)) {
          gain <- (table_per_cluster[k,2]-table_per_cluster[k+1,2])/table_per_cluster[k,2]
          if (gain < THRESHOLD[2]) {    # a very small one for this criterion
            break
          }
          print(paste("Gain", gain, "adding it to the DIC line"))
          gain_vector <- c(gain_vector, gain)
          x_elbow <- c(x_elbow, table_per_cluster[k+1,c("nclusters_pred")])
          y_elbow <- c(y_elbow, table_per_cluster[k+1,c("score")])
        }
        gfile <- file.path(outdir, paste0(flag, "_results_gain_", model, ".tsv"))
        gtable <- data.frame()
        gtable <- cbind(x_elbow, y_elbow, gain_vector)
        write.table(gtable, file=gfile, quote=FALSE, sep="\t", row.names=FALSE)

        print("x_elbow and y_elbow")
        print(x_elbow)
        print(y_elbow)
        print(length(x_elbow))
        # If there is only 1 point, that is the best cluster
        if (length(gain_vector) == 1) {
          best_nclusters <- x_elbow[1]
        }

        # SOmetimes there is more than one very large gain, remove the first
        if (length(gain_vector) >= 3) {
          if (gain_vector[2] >= 0.2 && gain_vector[3] >= 0.2) {
            x_elbow <- x_elbow[2:length(x_elbow)]
            y_elbow <- y_elbow[2:length(y_elbow)]
          }
        }
        print("x_elbow and y_elbow 2")
        print(x_elbow)
        print(y_elbow)

        # Now check that there is at least one gain that is >= THRESHOLD[1]
        # MA 10 May 2018: leaving THRESHOLD[1] below, but we may not need this to be different from the second threshold
        if (sum(gain_vector >= THRESHOLD[1]) >= 1) {
          # add one more point for elbow_finder
          x_elbow <- c(x_elbow, x_elbow[-1]+1)
          y_elbow <- c(y_elbow, y_elbow[-1])
          best_nclusters <- elbow_finder(x_elbow, y_elbow)
        } else {
          best_nclusters <- x_elbow[1]
        }
        best_row <- table_per_cluster[table_per_cluster$nclusters_pred == best_nclusters,]
        best_score <- table_per_cluster[table_per_cluster$nclusters_pred == best_nclusters, c("score")]
        print("Best row")
        print(best_row)
      }
    }

    ## Now plotting a graph that shows the DIC measures
    pdf(file.path(outdir, "DIC_selection.pdf"), height=7, width=9)
    x <- table_per_cluster[,c("nclusters_pred")]
    y <- table_per_cluster[,c("score")]

    # type='o' means it plots both points and lines overplotted
    matplot(x, y, lty=1, type='o', lwd=c(4), col=c(4), pch=19,
      ylab="DIC measure - proportion",
      xlab="Number of clusters",
      cex.axis=1.2,cex.lab=1.2)
    best_cluster <- best_row$nclusters_pred
    print(paste("Best cluster is", best_cluster))
    abline(v=best_cluster, col="red")
    if (criterion == "DIC_LINE_ELBOW") {
      abline(v=max(x_elbow), col="green")
    }
    axis(1, y)
    grid()
    dev.off()
  } else {
    print("Taking max")
    best_score <- max(table_all$score)
    best_row <- table_all[which.max(table_all$score),]
  }

  cpu_total <- sum(table_all$cpu_time)  # TOTAL CPU TIME
  memory <- best_row$memory  # MEMORY ONLY OF THE BEST RUN
  nclusters_pred <- best_row$nclusters_pred

  print (paste0('Selection criterion ', criterion))
  print (paste0('Total repeats ', ntotal))
  print (paste0('Number converged ', nconv))
  print (paste0('Best selection criterion ', best_score))
  best_cluster <- best_row$run
  print(paste("Best cluster", best_cluster))

  table_best <- data.frame(ntotal, nconv, best_score, best_cluster, cpu_total, memory, nclusters_pred)
  print ('Table of best results')
  print (table_best)

  return(list("best_row" = best_row, "table_best" = table_best))
}

compute_and_save_hd <- function (data_true, data_estimate, outfile, cpg_indicator_matrix = "") {
	# TODO: true data should be the original methylation file
	# first remove the columns that have NA in data_true
	# Or we can test by cell, but then when we sample we have to keep track of what we removed.

  data_true_orig <- data_true
	data_estimate <- data_estimate[,colSums(is.na(data_true))==0]
	data_true <- data_true[,colSums(is.na(data_true))==0]	

    ncells <- dim(data_true)[1]
    nloci <- dim(data_true)[2]
    hamming_distance_per_cell <- NULL
    
    print("number of cells")
    print(ncells)
    print(cpg_indicator_matrix)
    
    if (cpg_indicator_matrix == "") {
      print ("NO cpg_indicator_matri given, comparing all non-NA columns")
      for(i in 1:ncells) {
        hd <- sum(data_true[i,] != data_estimate[i,])/nloci
        hamming_distance_per_cell <- c(hamming_distance_per_cell,hd)
      }
      print("number of CpGs")
      print(nloci)
    } else {    # compare only the values that are 1 in the cpg_indicator_matrix
      print (paste0("Using cpg_indicator_matrix", cpg_indicator_matrix))
      cpg_ind <- read.csv(cpg_indicator_matrix, sep="\t", header=TRUE)
      cpg_ind <- as.matrix(cpg_ind[,-1])
      cpg_ind <- cpg_ind[,colSums(is.na(data_true_orig))==0]	
      colnames(data_true) <- colnames(cpg_ind)
      colnames(data_estimate) <- colnames(cpg_ind)
      for (i in 1:ncells) {
        mycolumns <- colnames(cpg_ind)[cpg_ind[i,]==1]
        hd <- sum(data_true[i,mycolumns] != data_estimate[i,mycolumns])/length(mycolumns)
        hamming_distance_per_cell <- c(hamming_distance_per_cell,hd)        
      }
    }

    # inferring the methylation profiles not from the G matrix, but from the clustering result



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
    return(hd_stats)
}


##### given a methylation file and a true clustering file, save the epigenotype file
## used in the subsampling for Smallwood
## Actually I shouldn't need this because I should use the whole full matrix
save.epigenotype <- function (methylation_file, true_Z_file, outfile) {
	library("matrixStats")
    true_Z <- read.csv(true_Z_file, sep="\t")
    meth <- read.csv(meth_file, sep="\t", header=TRUE)
    cluster1 <- true_Z[true_Z$epigenotype_id==1,]$cell_id
    cluster2 <- true_Z[true_Z$epigenotype_id==2,]$cell_id    
    meth1 <- meth[meth$cell_id %in% cluster1,]
    meth1$cell_id <- NULL
    epi1 <- c(1, colMedians(as.matrix(meth1), na.rm=TRUE))
    meth2 <- meth[meth$cell_id %in% cluster2,]
    meth2$cell_id <- NULL
    epi2 <- c(2, colMedians(as.matrix(meth2), na.rm=TRUE)) 
    epi <- t(data.frame(epi1,epi2))
    write.table (epi, file=outfile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

##====================================================================
#' Calculate hamming distance between true and estimated epigenotype/clusters
#'
#' @export
#'
#' @param outfile_est Path to the output file for the Bayesian estimates
#' @param true_epigenotype_file Path to true epigenotypes
#' @param true_Z_file Path to true cell membership
#' @param estimated_epigenotype_file Path to estimated epigenotypes
#' @param estimated_Z_file Path to estimated cell membership
#' @param methylation_file Path to the methylation file
#' @param regions_file Path to the region file
#'
#' @importFrom pheatmap pheatmap
#'
#' @examples
#' hamming.dist(outfile_est, true_epigenotype_file, true_Z_file, estimated_epigenotype_file, estimated_Z_file, methylation_file, regions_file)

hamming.dist <- function (outfile_est, true_epigenotype_file, true_Z_file, estimated_epigenotype_file, estimated_Z_file, methylation_file, regions_file, cpg_indicator_matrix="") {
  true_Z <- read.csv(true_Z_file, sep="\t")
  # TODO: for subsampling true data should be the original methylation file  ???
  ### Yes, add 
  # We don't usually have a header for this
  #if (!hd_from_full) {
    true_epi <- read.csv(true_epigenotype_file, sep="\t", header=FALSE)
    # the number of rows is the number of true clusters

    outfile_naive <- paste0(outfile_est, ".naive.tsv")
    outfile_est_corr <- paste0(outfile_est, ".corr.tsv")

    true_epi.tmp <- as.matrix(true_epi[,-1])
    data_true <- (t(sapply(true_Z[,2], function(x) return(true_epi.tmp[which(true_epi[,1]==x),])))) ### obtaining the true epigenotype for each cell
    rm(true_epi.tmp)
  #} else {
  #  # true_epigenotype_file is actually the full matrix
  #  data_true <- read.csv(true_epigenotype_file, sep="\t", header=TRUE)
  #}
  estimate_Z <- as.data.frame(read.csv(estimated_Z_file, sep='\t'))
  estimate_epi <- read.csv(estimated_epigenotype_file, sep='\t', header=TRUE)

  #########
  # compute the cellxCpG matrix of estimated CpGs
  estimate_epi.tmp <- as.matrix(estimate_epi[,-1])
  data_estimate <- (t(sapply(estimate_Z[,2], function(x) return(estimate_epi.tmp[which(estimate_epi[,1]==x),]))))
  rm(estimate_epi.tmp)

  # this will be the corrected version of the data estimate
  data_estimate_corr <- as.matrix(data_estimate)
  #########

  meth_data <- as.matrix(read.csv(methylation_file, sep='\t', header=TRUE, row.names=1))
  regions <- read.csv(regions_file, sep='\t', header=TRUE, row.names=1)
  num_regions <- dim(regions)[1]

  # impute the missing meth data the naive way
  # we used to traverse by CpG, but now we traverse by region and then by CpG

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
      variable_region <- .find_variable_region(r, cluster_set, estimate_epi, rstart, rend)
    }
    .impute_missing_data(meth_data, data_estimate_corr, estimate_Z, rstart, rend, variable_region)
  }
  
  print("1")
  estimates <- compute_and_save_hd(data_true, data_estimate, outfile_est, cpg_indicator_matrix=cpg_indicator_matrix)
  print("2")
  corrected_estimates <- compute_and_save_hd(data_true, data_estimate_corr, outfile_est_corr, cpg_indicator_matrix=cpg_indicator_matrix)
    print("3")
  naive_data <- compute_and_save_hd(data_true, meth_data, outfile_naive, cpg_indicator_matrix=cpg_indicator_matrix)
    print("4")

  return(list("estimates" = estimates, 'corrected_estimates' = corrected_estimates, 'naive_data' = naive_data))
}
