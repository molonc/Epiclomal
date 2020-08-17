
suppressMessages(library("argparse"))
suppressMessages(library("yaml"))
suppressMessages(library("REpiclomal"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

parser$add_argument("--input_dir", type="character", help="Input directory of epiclomal results")
parser$add_argument("--output_dir", type="character", help="Output directory of the evaluation results")
parser$add_argument("--model_name", type="character", help="A name for the model")
parser$add_argument("--methylation_file", type="character", default=NULL, help="Input data methylation file")
parser$add_argument("--regions_file", type="character", default=NULL, help="Input regions file, has to be given even for basic")
parser$add_argument("--true_clusters_file", type="character", default=NULL, help="File with the true clusters, if known")
parser$add_argument("--true_epigenotypes_file", type="character", default=NULL, help="File with the true epigenotypes, if known")
parser$add_argument("--cpg_indicator_matrix", type="character", default="", help="A file with a cell x cpg matrix ")
# this is used in the subsampling

# GAIN_THRESHOLD <- 0.05
args <- parser$parse_args()
print(args)

# args <- commandArgs(TRUE)

input <- args$input_dir
outdir <- args$output_dir
model <- args$model_name
meth_file <- args$methylation_file
regions_file <- args$regions_file
true_clusters_file <- args$true_clusters_file  # NULL if it is not known
true_epigenotypes_file <- args$true_epigenotypes_file  # NULL if it is not known
cpg_indicator_matrix <- args$cpg_indicator_matrix   # NULL it not given

######################################

run_eval <- function (input, flag, criterion, GAIN_THRESHOLD) {
  print(paste("criterion:", criterion, "GAIN_THRESHOLD:", GAIN_THRESHOLD))
  evaluation <- evaluate.epiclomal (input, outdir, model, flag, criterion, GAIN_THRESHOLD)  # essentially means the elbow curve can go all the way to the end
  best_row <- evaluation$best_row
  table_best <- evaluation$table_best
  best_cluster <- best_row$run
  epiMAPfile <- gsub("cluster_posteriors.tsv.gz", "genotype_MAP.tsv.gz", best_cluster)
  clMAPfile  <- gsub("cluster_posteriors.tsv.gz", "cluster_MAP.tsv.gz", best_cluster)
  # We add best_vmeasure later, only if there is a true_clusters_file

  outdir <- file.path(outdir, paste0(criterion, "_gainthr", GAIN_THRESHOLD))
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  if (!is.null(true_clusters_file))
  {
      # get the true number of clusters
      true_clusters_all <- read.table(true_clusters_file, header=TRUE)
      true_clusters <- true_clusters_all[["epigenotype_id"]]
      nclusters_true <- length(unique(true_clusters))
      print(paste0("Number of true clusters: ", nclusters_true))

      best_vmeasure <- best_row$all_vmeasure
      clone_prev_MAE <- best_row$clone_prev_MAE
      clone_prev_MSE <- best_row$clone_prev_MSE
      slsbulk_vmeasure <- best_row$slsbulk_vmeasure
      slsbulk_clone_prev_MAE <- best_row$slsbulk_clone_prev_MAE
      slsbulk_clone_prev_MSE <- best_row$slsbulk_clone_prev_MSE
      uncertainty <- best_row$uncertainty

      print(paste0('Vmeasure is ', best_vmeasure))
      table_best <- cbind (table_best, nclusters_true, best_vmeasure, clone_prev_MAE, clone_prev_MSE, slsbulk_vmeasure, slsbulk_clone_prev_MAE, slsbulk_clone_prev_MSE, uncertainty)
      # print ('Table with vmeasure')
      # print (table_best)

      if (!is.null(true_epigenotypes_file))
      {
          houtfile <- file.path(outdir, paste0(flag, "_hdist_bestrun_", model, ".tsv"))

          print ("Calling the hamming distance software")
          hamming.dist(houtfile, true_epigenotypes_file, true_clusters_file, epiMAPfile, clMAPfile, meth_file, regions_file, cpg_indicator_matrix)

          # get the mean value and put it in the table below
          hvalues <- read.table(houtfile, header=TRUE)
          hmean1 <- hvalues["mean"]
          colnames(hmean1)<-"hd_mean"
          
          hvalues <- read.table(paste0(houtfile, ".corr.tsv"), header=TRUE)
          hmean2 <- hvalues["mean"]
          colnames(hmean2)<-"hd_corr_mean"    
        
          hvalues <- read.table(paste0(houtfile, ".naive.tsv"), header=TRUE)
          hmean3 <- hvalues["mean"]
          colnames(hmean3)<-"hd_naive_mean"                 

          print(paste0('HD mean is ', hmean1))
          table_best <- cbind (table_best, hmean1, hmean2, hmean3)
          print('Table after hdist')
          print(table_best)
      }
  }

  sfile <- file.path(outdir, paste0(flag, "_results_bestrun_", model, ".tsv"))
  write.table(table_best, file=sfile, quote = FALSE, sep = "\t", row.names = FALSE)

  # Order will be by predicted
  print ("Calling the visualization software")
  print(clMAPfile)
  visualization(outdir=outdir,
    input_CpG_data_file=meth_file,
    input_regions=regions_file,
    inferred_clusters_file=clMAPfile,
    true_clusters_file=true_clusters_file,
    order_by_true=0,
    name=paste0(flag, "_bestrun_order_pred"))


  # if we know the true clusters, also print plots in which order is by true
  if (!is.null(true_clusters_file))
  {
      print ("Calling the visualization software the second time")
      visualization(outdir=outdir,
        input_CpG_data_file=meth_file,
        input_regions=regions_file,
        inferred_clusters_file=clMAPfile,
        true_clusters_file=true_clusters_file,
        order_by_true=1,
        name=paste0(flag, "_bestrun_order_true"))
  }
}
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



