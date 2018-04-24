# plotting Vmeasure, Hamming Distance, number of clusters, cluster prevalence mean absolute error and mean squared error

#.libPaths(c("/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.2", 
#            "/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.2", "/clusterapp/software/linux-x86_64-centos6/R-3.2.3/lib64/R/library","/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.3"))



suppressMessages(library("argparse"))
library(ggplot2)
library(gridExtra)
library(stringr)
library(plyr)
# library(ggpubr,lib.loc = "/clusterapp/clusterhome/csouza/R/x86_64-pc-linux-gnu-library/3.3")


# function to get the script path
getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}



scriptPath <- getScriptPath()
source(paste0(scriptPath, "/plot_functions.R"))
#library(ggpubr)

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

parser$add_argument("--output_dir", type="character", default="output", help="Entire or just the beginning of the output directory file ")

parser$add_argument("--criterion", type="character", help="The selection criterion for the best run: DIC_measure_gainthr0.05 DIC_measure_gainthr0.1")

datasets <- c("InHouse",
              "Smallwood2014",
              "Hou2016",
              "Luo2017",
              "Farlik2016")
datapaths <- c("../EPI-112_inhouse_data/FINAL_RESULTS",
               "../EPI-70_Smallwood2014/FINAL_RESULTS_IQR",
               "../EPI-105_scTrio/FINAL_RESULTS",
               "../EPI-106_Luo2017/FINAL_RESULTS",
               "../EPI-89_Farlik2016_all_union_IQR/FINAL_RESULTS")
simplepaths <- c("../EPI-112_inhouse_data/OUTPUT_epiclomal_INHOUSE/RUN/epiclomal_INHOUSE_",
                 "../EPI-70_Smallwood2014/OUTPUT_epiclomal_Smallwood2014_IQR/RUN/epiclomal_Smallwood2014_",
                 "../EPI-105_scTrio/OUTPUT_epiclomal_scTrio/RUN/epiclomal_scTrio_",
                 "../EPI-106_Luo2017/OUTPUT_epiclomal_Luo2017_genebodies_500_clean_random_cells/RUN/epiclomal_Luo2017_genebodies_500_clean_random_cells_",
                 "../EPI-89_Farlik2016_all_union_IQR/OUTPUT_epiclomal_Farlik2016_all_union/RUN/epiclomal_Farlik2016_all_union_")

# each replicate file should be in inputs, for example inputs/Smallwood2014_replicates.txt

args <- parser$parse_args() 

print(args)

outdir <- args$output_dir
dir.create(outdir, showWarnings = FALSE)

criterion <- args$criterion
dir.create(paste0(outdir,"/",criterion), showWarnings = FALSE)
outdir <- paste0(outdir,"/",criterion)


##################
### box plots and line plots ####
##################


##################
### plots clone_prev_MAE ####
##################
print ("Plots for clone_prev_MAE")
model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
ourcolors <- c("red", "blue", "green", "purple", "cyan")
label <- c("Epiclomal","EuclideanClust","DensityCut","HammingClust","PearsonClust")
plot_data(model, label, ourcolors, criterion, "clone_prev_MAE")

##################
### plots V-measure ####
##################
### V-measure
print ("Plots for V-measure")
model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
ourcolors <- c("red", "blue", "green", "purple", "cyan")
label <- c("Epiclomal","EuclideanClust","DensityCut","HammingClust","PearsonClust")
plot_data (model, label, ourcolors, criterion, "Vmeasure")


##################
### plots nclusters ####
##################
print ("Plots for nclusters")
model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
ourcolors <- c("red", "blue", "green", "purple", "cyan")
label <- c("Epiclomal","EuclideanClust","DensityCut","HammingClust","PearsonClust")
plot_data (model, label, ourcolors, criterion, "nclusters")

