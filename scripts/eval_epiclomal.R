
suppressMessages(library("argparse"))

RSCRIPT <- Sys.getenv("RSCRIPT")

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

parser$add_argument("--input_dir", type="character", help="Input directory of epiclomal results")
parser$add_argument("--output_dir", type="character", help="Output directory of the evaluation results")
parser$add_argument("--model_name", type="character", help="A name for the model")
parser$add_argument("--hdist_software", type="character", default=NULL, help="Full path of software that calculates hamming distance")
parser$add_argument("--visualization_software", type="character", default=NULL, help="Full path of software for visualization")
parser$add_argument("--methylation_file", type="character", default=NULL, help="Input data methylation file")
parser$add_argument("--regions_file", type="character", default=NULL, help="Input regions file, has to be given even for basic")
parser$add_argument("--true_clusters_file", type="character", default=NULL, help="File with the true clusters, if known")
parser$add_argument("--true_epigenotypes_file", type="character", default=NULL, help="File with the true epigenotypes, if known")
# GAIN_THRESHOLD <- 0.05
args <- parser$parse_args()
print(args)

# args <- commandArgs(TRUE)

input <- args$input_dir
output <- args$output_dir
model <- args$model_name
hdist_software <- args$hdist_software
visualization_software <- args$visualization_software
meth_file <- args$methylation_file
regions_file <- args$regions_file
true_clusters_file <- args$true_clusters_file  # NULL if it is not known
true_epigenotypes_file <- args$true_epigenotypes_file  # NULL if it is not known

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

######################################

run_eval <- function (input, flag, criterion, GAIN_THRESHOLD)  {
    # criterion can be "lower_bound" or "log_posterior"
    output <- paste0(output,"/",criterion,"_gainthr", GAIN_THRESHOLD)
    dir.create(output, showWarnings=FALSE, recursive=TRUE)

    # input can be actually a list of directories. Then look through all of them and compute the measure
    directories <- Sys.glob(input)
    files <- NULL
    # NOTE: This may be very slow for many files
    for (dir in directories) {
        files <- c(files, list.files(dir, recursive = TRUE, pattern = "params.yaml", full.names = TRUE))
    }

    lines <- lapply(files, readLines)

    # print(lines)

    run <- NULL
    # NOTE: This may be very slow for many files
    for (dir in directories) {
        run <- c(run, list.files(dir, recursive = TRUE, pattern = "cluster_posteriors.tsv.gz", full.names = TRUE))
    }
    #print (run)
    # dir <- paste0(getwd(), "/", input, "/", run)
    converged <- sapply(sapply(lines, grep, pattern = "converged: true", value = TRUE), length)
    criterion <- paste0(criterion, ": ")
    if (criterion == "DIC_measure: " || criterion == "DIC_LINE_ELBOW: ")  {
        measure <- "DIC_measure: "
    } else {
        measure <- criterion
    }

    score <- as.numeric(sub(measure, "", sapply(lines, grep, pattern = measure, value = TRUE)))
    # elbo <- as.numeric(sub(paste0(criterion,": "), "", sapply(lines, grep, pattern = criterion, value = TRUE)))
    # Now elbo variable can be the elbo (lower_bound) or the log_posterior (unnormalized)
    cpu_time <- as.numeric(sub("CPU_time_seconds: ", "", sapply(lines, grep, pattern = "CPU_time_seconds", value = TRUE)))
    memory <- as.numeric(sub("Max_memory_MB: ", "", sapply(lines, grep, pattern = "Max_memory_MB", value = TRUE)))
    all_vmeasure <- as.numeric(sub("Vmeasure: ", "", sapply(lines, grep, pattern = "Vmeasure", value = TRUE)))
    nclusters_pred <- as.numeric(sub("nclusters: ", "", sapply(lines, grep, pattern = "nclusters", value = TRUE)))
    # for clone_prev_MAE, I may also have slsbulk_clone_prev_MAE
    clone_prev_MAE <- as.numeric(sub("clone_prev_MAE: ", "", sapply(lines, grep, pattern = "^clone_prev_MAE", value = TRUE)))
    clone_prev_MSE <- as.numeric(sub("clone_prev_MSE: ", "", sapply(lines, grep, pattern = "^clone_prev_MSE", value = TRUE)))
    slsbulk_vmeasure <- as.numeric(sub("slsbulk_vmeasure: ", "", sapply(lines, grep, pattern = "slsbulk_vmeasure", value = TRUE)))
    slsbulk_clone_prev_MAE <- as.numeric(sub("slsbulk_clone_prev_MAE: ", "", sapply(lines, grep, pattern = "slsbulk_clone_prev_MAE", value = TRUE)))
    slsbulk_clone_prev_MSE <- as.numeric(sub("slsbulk_clone_prev_MSE: ", "", sapply(lines, grep, pattern = "slsbulk_clone_prev_MSE", value = TRUE)))
    uncertainty <- as.numeric(sub("uncertainty_true_positive_rate: ", "", sapply(lines, grep, pattern = "uncertainty_true_positive_rate", value = TRUE)))

    table_all <- data.frame(converged, score, run, cpu_time, memory, nclusters_pred, all_vmeasure, clone_prev_MAE, clone_prev_MSE, slsbulk_vmeasure, slsbulk_clone_prev_MAE, slsbulk_clone_prev_MSE, uncertainty)

    dfile <- paste0(output, "/", flag, "_results_allruns_", model, ".tsv")

    write.table(table_all, file = dfile, quote = FALSE, sep = "\t", row.names = FALSE)


    sfile <- paste0(output, "/", flag, "_results_bestrun_", model, ".tsv")

    ntotal <- nrow(table_all)
    nconv <- sum(table_all[,1])
    if (criterion == "DIC_measure: " || criterion == "DIC_LINE_ELBOW: ")  {

        # print('Taking min')
        # best_score <- min(table_all[,2])
        # bestrow <- which.min(table_all[,2])

        table_per_cluster = data.frame()
        # find the unique number of clusters
        for ( k in sort(unique(table_all[,6])))   {
            # print (k)
            rows <- table_all[which(table_all[,6]==k),]
            # print(rows)
            minofk_score <- min(rows[,2])
            minofk_row <- rows[which.min(rows[,2]),]
            # print(minofk_score)
            # print(minofk_row)
            table_per_cluster = rbind(table_per_cluster, minofk_row)
        }
        # Now measure the relative gain
        print("Table per cluster")
        print(table_per_cluster)
        dfile <- paste0(output, "/", flag, "_results_perclusterruns_", model, ".tsv")
        write.table(table_per_cluster, file = dfile, quote = FALSE, sep = "\t", row.names = FALSE)

        # also plot the v-measure versus the number of predicted clusters to see if we get better V-measure when we choose a different number of clusters
        pdf(paste0(output,"/Vmeasure_vs_nclusters.pdf"),height=7,width=9)
        x <- table_per_cluster[,c("nclusters_pred")]
        y <- table_per_cluster[,c("all_vmeasure")]

        # type='o' means it plots both points and lines overplotted
        matplot(x,y,lty=1,type='o',lwd=c(4),col=c(4), pch=19,
            ylab="V-measure for run with best score",
            xlab="Number of clusters",
            cex.axis=1.2,cex.lab=1.2)
        dev.off()

        print(paste0("Number of rows in table per cluster is ", nrow(table_per_cluster)))

        if (nrow(table_per_cluster) == 1) {
            best_score <- table_per_cluster[1,2]
            bestrow <- table_per_cluster[1,]
        } else {
            if (criterion == "DIC_measure: ") {
                for ( k in 1:(nrow(table_per_cluster)-1)) {
                    gain <- (table_per_cluster[k,2]-table_per_cluster[k+1,2])/table_per_cluster[k,2]
                    print(paste0("Gain ", gain))
                    if (gain < GAIN_THRESHOLD) {
                        best_score <- table_per_cluster[k,2]
                        bestrow <- table_per_cluster[k,]
                        break
                    } else {
                        best_score <- table_per_cluster[k+1,2]
                        bestrow <- table_per_cluster[k+1,]
                    }
                    # TODO: If nothing was found, it means we have to run with more clusters
                }
            }
            if (criterion == "DIC_LINE_ELBOW: ")  {
                print ("Doing DIC LINE ELBOW")
                print (table_per_cluster)
                THRESHOLD <- as.numeric(unlist(strsplit(GAIN_THRESHOLD, split ="_")))
                # first add point 1
                x_elbow <- c(table_per_cluster[1,c("nclusters_pred")])
                y_elbow <- c(table_per_cluster[1,c("score")])
                gain_vector <- c(0)
                for ( k in 1:(nrow(table_per_cluster)-1)) {
                    gain <- (table_per_cluster[k,2]-table_per_cluster[k+1,2])/table_per_cluster[k,2]
                    if (gain < THRESHOLD[2]) {    # a very small one for this criterion
                        break
                    }
                    print(paste0("Gain ", gain, " adding it to the DIC LINE"))
                    gain_vector <- c(gain_vector, gain)
                    x_elbow <- c(x_elbow, table_per_cluster[k+1,c("nclusters_pred")])
                    y_elbow <- c(y_elbow, table_per_cluster[k+1,c("score")])
                }
                gfile <- paste0(output, "/", flag, "_results_gain_", model, ".tsv")
                gtable <- data.frame()
                gtable <- cbind(x_elbow, y_elbow, gain_vector)
                write.table(gtable, file = gfile, quote = FALSE, sep = "\t", row.names = FALSE)

                print ("x_elbow and y_elbow")
                print (x_elbow)
                print (y_elbow)
                print (length(x_elbow))
                # If there is only 1 point, that is the best cluster
                if (length(x_elbow) == 1) {
                    best_nclusters <- x_elbow[1]
                }

                # SOmetimes there is more than one very large gain, remove the first
                if (length(gain_vector) >= 3) {
                    if (gain_vector[2] >= 0.2 && gain_vector[3] >= 0.2) {
                        # eliminate the first one
                        x_elbow <- x_elbow[2:length(x_elbow)]
                        y_elbow <- y_elbow[2:length(y_elbow)]
                    }
                }
                print ("x_elbow and y_elbow 2")
                print (x_elbow)
                print (y_elbow)
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
                bestrow <- table_per_cluster[table_per_cluster$nclusters_pred == best_nclusters,]
                best_score <- table_per_cluster[table_per_cluster$nclusters_pred == best_nclusters,c("score")]
                print("Best row")
                print (bestrow)
            }
        }
        ## Now plotting a graph that shows the DIC measures

        pdf(paste0(output,"/DIC_selection.pdf"),height=7,width=9)
        x <- table_per_cluster[,c("nclusters_pred")]
        # y <- table_per_cluster[,c("score")/as.numeric(table_per_cluster[1,c("score")])]
        y <- table_per_cluster[,c("score")]

        # type='o' means it plots both points and lines overplotted
        matplot(x,y,lty=1,type='o',lwd=c(4),col=c(4), pch=19,
            ylab="DIC measure - proportion",
            xlab="Number of clusters",
            #xaxt="n",
            cex.axis=1.2,cex.lab=1.2)
        #points(x, cex = 1.5, col = "dark red")
        bestcluster <- bestrow[,6]
        print(paste0('Best cluster is ', bestcluster))
        abline(v=bestcluster, col="red")
        if (criterion == "DIC_LINE_ELBOW: ") {
            abline(v=max(x_elbow), col="green")
        }
        axis(1, y)
        grid()

        dev.off()


    } else {
        print('Taking max')
        best_score <- max(table_all[,2])
        bestrow <- table_all[which.max(table_all[,2]),]
    }
    cpu_total <- sum(table_all[,4])     # TOTAL CPU TIME
    memory <- bestrow[,5]      # MEMORY ONLY OF THE BEST RUN
    nclusters_pred <- bestrow[,6]

    print (paste0('Selection criterion ', criterion))
    print (paste0('Total repeats ', ntotal))
    print (paste0('Number converged ', nconv))
    print (paste0('Best selection criterion ', best_score))
    bestcluster <- bestrow[,3]
    print (paste0('Best cluster ', bestcluster))

    epiMAPfile <- gsub("cluster_posteriors.tsv.gz", "genotype_MAP.tsv.gz", bestcluster)
    clMAPfile  <- gsub("cluster_posteriors.tsv.gz", "cluster_MAP.tsv.gz", bestcluster)
    # We add best_vmeasure later, only if there is a true_clusters_file
    table_best <- data.frame(ntotal, nconv, best_score, bestcluster, cpu_total, memory, nclusters_pred)
    print ('Table of best results')
    print (table_best)

    # Now call the visualization software
    visline <- paste0("--out_directory=", output,
                    " --methylation_file=", meth_file,
                    " --regions_file=", regions_file,
                    " --inferred_clusters_file=", clMAPfile)

    if (!is.null(true_clusters_file))
    {
        # get the true number of clusters
        true_clusters_all <- read.table(true_clusters_file, header=TRUE)
        true_clusters <- true_clusters_all[["epigenotype_id"]]
        nclusters_true <- length(unique(true_clusters))
        print(paste0("Number of true clusters: ", nclusters_true))

        best_vmeasure <- bestrow[,7]
        clone_prev_MAE <- bestrow[,8]
        clone_prev_MSE <- bestrow[,9]
        slsbulk_vmeasure <- bestrow[,10]
        slsbulk_clone_prev_MAE <- bestrow[,11]
        slsbulk_clone_prev_MSE <- bestrow[,12]
        uncertainty <- bestrow[,13]

        print(paste0('Vmeasure is ', best_vmeasure))
        table_best <- cbind (table_best, nclusters_true, best_vmeasure, clone_prev_MAE, clone_prev_MSE, slsbulk_vmeasure, slsbulk_clone_prev_MAE, slsbulk_clone_prev_MSE, uncertainty)
        # print ('Table with vmeasure')
        # print (table_best)

        visline <- paste0(visline, " --true_clusters_file=", true_clusters_file)

        if (!is.null(true_epigenotypes_file))
        {
            houtfile <- paste0(output, "/", flag, "_hdist_bestrun_", model, ".tsv")

            # setting up the arguments of hamming_distance.R
            hargs <- NULL
            hargs$output_file <- houtfile
            hargs$true_epigenotype_file <- true_epigenotypes_file
            hargs$true_membership_file <- true_clusters_file
            hargs$estimated_epigenotype_file <- epiMAPfile
            hargs$estimated_membership_file <- clMAPfile
            args <- hargs
            hline <- paste0("--output_file=", houtfile, " --true_epigenotype_file=", true_epigenotypes_file,
                    " --true_membership_file=", true_clusters_file, " --estimated_epigenotype_file=", epiMAPfile,
                    " --estimated_membership_file=", clMAPfile, " --methylation_file=", meth_file,
                    " --regions_file=", regions_file)
            print ("Calling the hamming distance software")
            command <- paste (RSCRIPT, hdist_software, hline)
            print(command)
            system(command)
            # source (hdist_software)

            # get the mean value and put it in the table below
            hvalues <- read.table(houtfile, header=TRUE)
            hmean <- hvalues["mean"]
            colnames(hmean)<-"hd_mean"

            print(paste0('HD mean is ', hmean))
            table_best <- cbind (table_best, hmean)
            print('Table after hdist')
            print(table_best)
        }
    }

    write.table(table_best, file=sfile, quote = FALSE, sep = "\t", row.names = FALSE)

    # Order will be by predicted
    print ("Calling the visualization software")
    command <- paste0 (RSCRIPT, " ", visualization_software, " ", visline, " --order_by_true=0 --name=", flag, "_bestrun_order_pred")
    print(command)
    system(command)


    # if we know the true clusters, also print plots in which order is by true
    if (!is.null(true_clusters_file))
    {
        print ("Calling the visualization software the second time")
        command <- paste0 (RSCRIPT, " ", visualization_software, " ", visline, " --order_by_true=1 --name=", flag, "_bestrun_order_true")
        #command <- paste0 ("Rscript ", visualization_software, " ", visline, " --order_by_true=1 --name=InHouse")
        print(command)
        system(command)
    }
}   # end function run_eval

######################################



# run_eval (input, "all", "lower_bound")
# run_eval (input, "all", "log_posterior_allK")
# run_eval (input, "all", "log_posterior_clusterK")
# run_eval (input, "all", "log_likelihood")
# run_eval (input, "all", "DIC_measure", 0.03)
# run_eval (input, "all", "DIC_measure", 0.1)

run_eval (input, "all", "DIC_LINE_ELBOW", "0.05_-100")  # essentially means the elbow curve can go all the way to the end
run_eval (input, "all", "DIC_LINE_ELBOW", "0.05_0.02")
run_eval (input, "all", "DIC_LINE_ELBOW", "0.02_0.02")
run_eval (input, "all", "DIC_measure", 0.05)

#run_eval (input, "all", "DIC_LINE_ELBOW", "0.02_0.001")
#
# run_eval (input, "all", "DIC_LINE_ELBOW", "0.1_0.02")

#input <- paste0(input, "/{0..9}/")
#print(input)
#run_eval (input, "init_hdist")



