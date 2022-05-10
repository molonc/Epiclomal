# plotting Vmeasure, Hamming Distance, number of clusters, cluster prevalence mean absolute error and mean squared error

#.libPaths(c("/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.2", 
#            "/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.2", "/clusterapp/software/linux-x86_64-centos6/R-3.2.3/lib64/R/library","/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.3"))



suppressMessages(library("argparse"))
library(stringr)
# library(ggpubr,lib.loc = "/clusterapp/clusterhome/csouza/R/x86_64-pc-linux-gnu-library/3.3")


# function to get the script path
# getScriptPath <- function(){
#     cmd.args <- commandArgs()
#     m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
#     script.dir <- dirname(regmatches(cmd.args, m))
#     if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
#     if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
#     return(script.dir)
# }
# 
# scriptPath <- getScriptPath()

#source(paste0(scriptPath, "/plot_functions.R"))
#library(ggpubr)

source("/scratch/shahlab_tmp/mandronescu/Epiclomal/Epiclomal/scripts/plot_functions.R")

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

parser$add_argument("--output_dir", type="character", default="output", help="Entire or just the beginning of the output directory file ")

parser$add_argument("--criterion1", type="character", help="The non-very granular selection criterion")
parser$add_argument("--criterion2", type="character", help="The granular selection criterion, to be used for Luo and Farlik")

datasets <- c(
#              "InHouse",
              "Smallwood2014",
              "Hou2016",
#              "Luo2017",
              "Farlik2016")
ntrue_clusters <- c(3,2,2,21,6)              
              
datapaths <- c(
#               "/shahlab/mandronescu/EPI-112_inhouse_data/FINAL_RESULTS",
               "/shahlab/mandronescu/EPI-70_Smallwood2014/FINAL_RESULTS",
               "/shahlab/mandronescu/EPI-105_scTrio/FINAL_RESULTS",
#               "/shahlab/mandronescu/EPI-106_Luo2017/FINAL_RESULTS_genebodies_all_clean_cells_MAXK30",
               "/shahlab/mandronescu/EPI-89_Farlik2016/FINAL_RESULTS_6clusters_Farlik_clustering")
simplepaths <- c(
#                 "/shahlab/mandronescu/EPI-112_inhouse_data/OUTPUT_epiclomal_INHOUSE/RUN/epiclomal_INHOUSE_",
                 "/shahlab/mandronescu/EPI-70_Smallwood2014/OUTPUT_epiclomal_Smallwood2014/RUN/epiclomal_Smallwood2014_",
                 "/shahlab/mandronescu/EPI-105_scTrio/OUTPUT_epiclomal_scTrio/RUN/epiclomal_scTrio_",
#                 "/shahlab/mandronescu/EPI-106_Luo2017/OUTPUT_epiclomal_Luo2017_genebodies_all_clean_cells_MAXK30/RUN/epiclomal_Luo2017_genebodies_all_clean_cells_MAXK30_",
                 "/shahlab/mandronescu/EPI-89_Farlik2016/OUTPUT_epiclomal_Farlik2016_6clusters_Farlik_clustering/RUN/epiclomal_Farlik2016_6clusters_Farlik_clustering_")

# each replicate file should be in inputs, for example inputs/Smallwood2014_replicates.txt

args <- parser$parse_args() 

print(args)

outdir <- args$output_dir
dir.create(outdir, showWarnings = FALSE)

criterion1 <- args$criterion1
criterion2 <- args$criterion2
dir.create(paste0(outdir,"/",criterion1,"__",criterion2), showWarnings = FALSE)
outdir <- paste0(outdir,"/",criterion1,"__",criterion2)

all_regions_criterion <- "0_1_0.01"


#################################
### box plots and line plots ####
#################################


collect_data <- function(model, criterion1, criterion2, measure_name){
  # measure_name can be HD, Vmeasure, nclusters, cp_error
  
  if (measure_name == "HD") {
    column <- "mean"
  } else if (measure_name == "Vmeasure") {
    column <- "best_vmeasure"
  } else if (measure_name == "nclusters") {
    column <- "nclusters_pred"
  } else if (measure_name == "nclusters2") {
    column <- "nclusters_pred"
  } else if (measure_name == "clone_prev_MAE") {
    column <- "clone_prev_MAE"
  } else if (measure_name == "clone_prev_MSE") {
    column <- "clone_prev_MSE"
  }    
  
  variable <- datasets
  # variable is the value of the changed variable, for example if we are varying misspb, variable is 0.5, 0.6, 0.7, 0.8, 0.9, 0.95
  
  savedfile <- paste0(outdir,"/data_",measure_name,"_",criterion1,"__",criterion2,".Rda")
  #print(savedfile) 
  
 if (file.exists(savedfile)) {
   print("File already exists, loading it")
   load(savedfile)
 } else {
    print("File doesn't exist, creating it")
    method <- NULL
    VAR <- NULL
    measure <- NULL

    replicate <- NULL
    imputed <- NULL
    
    crash <- NULL
    
    for(m in 1:length(model)){
      for(j in 1:length(variable)){
        replicate_file <- paste0("inputs/", variable[j], "_replicates.txt")
        replicates <- read.table (replicate_file, header=FALSE, sep="\t")
        
        #print(replicates)
  
        number_replicates <- nrow(replicates)
        print(paste0("Model ", model[m], " data set ", variable[j], " number of replicates ", number_replicates))
        
        for(i in 1:number_replicates){
        if (model[m] == "HammingClust" || model[m] == "DensityCut" || model[m] == "PearsonClust" || model[m] == "EuclideanClust") {
          #if (model[m] == "PBALclust" || model[m] == "densitycut" || model[m] == "Pearsonclust" || model[m] == "Hclust") {
            if (length(grep(all_regions_criterion, replicates[i,])) > 0) {
                results_file <- paste0(simplepaths[j],replicates[i,],"_epiclomal_real_task1/outputs/simple_hclust/results_", model[m], ".txt")
            }                
            else {
                results_file <- paste0(simplepaths[j],replicates[i,],"_epiclomal_real/outputs/simple_hclust/results_", model[m], ".txt")
            }                
          } else {   
            # Use the second (granular) criterion for Luo and Farlik
            if (datasets[j] == "Luo2017" | datasets[j] == "Farlik2016") {  
                rep <- replicates[i,]
                rep <- gsub("_K10","",rep)
                rep <- gsub("_K30","",rep)
                results_file <- paste0(datapaths[j],"/",rep,"_",model[m],"/",criterion2,"/all_results_bestrun_",model[m],".tsv")
            } else {    
                rep <- gsub("_K10","",replicates[i,])
                results_file <- paste0(datapaths[j],"/",rep,"_",model[m],"/",criterion1,"/all_results_bestrun_",model[m],".tsv")
            }
          }    
          print (paste0('Results file is ', results_file))                    
          
          t <- try(read.table(file=results_file,sep="\t",header=TRUE))   
          if("try-error" %in% class(t)) { 
                ### could have an alternativeFunction() here
                print("can't find file, trying with imputation")
                results_file <- gsub("simple_hclust", "simple_hclust_imputed", results_file)
                print (paste0('Imputed Results file is ', results_file))
                
                # Adding a star to label to show that it is using imputed values
                # replicates[i,] <- paste0(replicates[i,]," *")
                
                imputed <- c(imputed, " *")
                
                t <- try(read.table(file=results_file,sep="\t",header=TRUE))   
                if("try-error" %in% class(t)) {     # neither non-imputed nor imputed exists
                    print("can't find impute file")
                    crash <- c(crash,0)            
                    measure <- c(measure,NA)            
                    VAR <- c(VAR,variable[j])
                    method <- c(method,model[m]) 
                    replicate <- c(replicate,as.character(replicates[i,])) 
                    print(method)
                    print(replicate)                
                } else {  # imputed variant exists
                    print ("Found imputed file")
                    f <- read.table(file=results_file,sep="\t",header=TRUE)
                    this_measure <- f[,column]
                    # correct for the number of cluster differential
                    if (measure_name == "nclusters2") {
                        this_measure <- f[,column] - ntrue_clusters[j]
                    }                                        
                    crash <- c(crash,1)            
                    measure <- c(measure, this_measure)
                    VAR <- c(VAR,variable[j])
                    method <- c(method,model[m]) 
                    replicate <- c(replicate,as.character(replicates[i,])) 
                    print(method)
                    print(replicate)   
                }                            
          } else {
                f <- read.table(file=results_file,sep="\t",header=TRUE)
                this_measure <- f[,column]
                # correct for the number of cluster differential
                if (measure_name == "nclusters2") {
                    this_measure <- f[,column] - ntrue_clusters[j]
                }                   
                crash <- c(crash,1)            
                measure <- c(measure, this_measure)
                VAR <- c(VAR,variable[j])
                method <- c(method,model[m]) 
                replicate <- c(replicate,as.character(replicates[i,]))            
                imputed <- c(imputed, "")
                print(method)
                print(replicate)          
          }      
        }
      }
    }
    
    # make the measure be at least some small value so we can see it
    measure[measure<0.02] <- 0.02
    
    big_df <- cbind(as.data.frame(measure),as.data.frame(crash),VAR,method,replicate, imputed)
    colnames(big_df) <- c("Measure","crash","VAR","method","replicate","imputed")

    big_df$replicate <- as.character(big_df$replicate)
    
    str(big_df)
    
    big_df$VAR <- factor(big_df$VAR,levels=variable)
    
    print(str(big_df))
       
    # Now saving the data frame
    save(big_df, crash, file=savedfile)
  }  # end make the data files  
  return (list("big_df"=big_df, "crash"=crash))
}  



xlsfile <- paste0(outdir,"/SourceData.xlsx")
if (file.exists(xlsfile)) {
    file.remove(xlsfile)
} 

#############################
### plots clone_prev_MAE ####
#############################

print ("Plots for clone_prev_MAE")
#model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
model <- c("region", "EuclideanClust", "DensityCut", "HammingClust", "PearsonClust")
mylist <- collect_data(model, criterion1, criterion2, "clone_prev_MAE")
# plot_data(mylist$big_df, mylist$crash, model, "clone_prev_MAE",add_points=TRUE)
plot_data_barplots(mylist$big_df, mylist$crash, model,"clone_prev_MAE")

########################
### plots V-measure ####
########################

print ("Plots for V-measure")
#model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
model <- c("region", "EuclideanClust", "DensityCut", "HammingClust", "PearsonClust")
mylist <- collect_data (model, criterion1, criterion2, "Vmeasure")
# plot_data(mylist$big_df, mylist$crash, model, "Vmeasure",add_points=TRUE)
plot_data_barplots(mylist$big_df, mylist$crash, model,"Vmeasure")

########################
### plots nclusters ####
########################
print ("Plots for nclusters")
#model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
model <- c("region", "EuclideanClust", "DensityCut", "HammingClust", "PearsonClust")
# ourcolors <- c("red", "blue", "green", "purple", "cyan")
mylist <- collect_data (model, criterion1, criterion2, "nclusters")
#print(mylist$big_df)
# plot_data(mylist$big_df, mylist$crash, model, "nclusters",add_points=TRUE)
plot_data_barplots(mylist$big_df, mylist$crash, model,"nclusters")

########################
### plots nclusters ####
########################
print ("Plots for nclusters2")
#model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
model <- c("region", "EuclideanClust", "DensityCut", "HammingClust", "PearsonClust")
# ourcolors <- c("red", "blue", "green", "purple", "cyan")
mylist <- collect_data (model, criterion1, criterion2, "nclusters2")
#print(mylist$big_df)
# plot_data(mylist$big_df, mylist$crash, model, "nclusters",add_points=TRUE)
plot_data_barplots(mylist$big_df, mylist$crash, model,"nclusters2")



