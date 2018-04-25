# plotting Vmeasure, Hamming Distance, number of clusters, cluster prevalence mean absolute error and mean squared error

#.libPaths(c("/home/mandronescu/R/x86_64-pc-linux-gnu-library/3.2", 
#            "/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.2", "/clusterapp/software/linux-x86_64-centos6/R-3.2.3/lib64/R/library","/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos6/3.3"))



suppressMessages(library("argparse"))
library(stringr)
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



collect_data <- function(model, label, ourcolors, criterion, measure_name){
  # measure_name can be HD, Vmeasure, nclusters, cp_error
  
  if (measure_name == "HD") {
    column <- "mean"
  } else if (measure_name == "Vmeasure") {
    column <- "best_vmeasure"
  } else if (measure_name == "nclusters") {
    column <- "nclusters_pred"
  } else if (measure_name == "clone_prev_MAE") {
    column <- "clone_prev_MAE"
  } else if (measure_name == "clone_prev_MSE") {
    column <- "clone_prev_MSE"
  }    
  
  variable <- datasets
  # variable is the value of the changed variable, for example if we are varying misspb, variable is 0.5, 0.6, 0.7, 0.8, 0.9, 0.95
  
  savedfile <- paste0(outdir,"/data_",measure_name,"_",criterion,".Rda")
  if (file.exists(savedfile)) {
    print("File already exists, loading it")
    load(savedfile)
  } else {
    print("File doesn't exist, creating it")
    method <- NULL
    VAR <- NULL
    measure <- NULL
    
    crash <- NULL
    
    for(m in 1:length(model)){
      for(j in 1:length(variable)){
        replicate_file <- paste0("inputs/", variable[j], "_replicates.txt")
        replicates <- read.table (replicate_file, header=FALSE, sep="\t")
        number_replicates <- nrow(replicates)
        print(paste0("Model ", model[m], " data set ", variable[j], " number of replicates ", number_replicates))
        for(i in 1:number_replicates){
          if (model[m] == "PBALclust" || model[m] == "densitycut" || model[m] == "Pearsonclust" || model[m] == "Hclust") {
            results_file <- paste0(simplepaths[j],replicates[i,],"_K1_epiclomal_real/outputs/simple_hclust/results_", model[m], ".txt")                        
          } else {     
            results_file <- paste0(datapaths[j],"/",replicates[i,],"_",model[m],"/",criterion,"/all_results_bestrun_",model[m],".tsv")
          }    
          print (paste0('Results file is ', results_file))
          t <- try(read.table(file=results_file,sep="\t",header=TRUE))   
          if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
            print("can't find file")
            crash <- c(crash,0)            
            measure <- c(measure,NA)            
            VAR <- c(VAR,variable[j])
            method <- c(method,model[m]) 
            
          } else {
            f <- read.table(file=results_file,sep="\t",header=TRUE)
            crash <- c(crash,1)            
            measure <- c(measure,f[,column])
            VAR <- c(VAR,variable[j])
            method <- c(method,model[m]) 
          }      
        }
      }
    }
    
    
    big_df <- cbind(as.data.frame(measure),as.data.frame(crash),VAR,method)
    colnames(big_df) <- c("Measure","crash","VAR","method")
    str(big_df)
    
    ### changing variable names
    
    for (i in 1:length(model)){
      
      big_df$method <- sub(pattern=model[i],x=big_df$method,replacement=label[i])
      
    }
    
    big_df$method <- factor(big_df$method,levels=label)
    
    big_df$VAR <- factor(big_df$VAR,levels=variable)
    
    print("Big DF")
    
    print(big_df)
    
    print(str(big_df))
    
    
    # Now saving the data frame
    save(big_df, crash, file=savedfile)
  }  # end make the data files  
  return (list(big_df, crash))
}  

##################
### plots clone_prev_MAE ####
##################
print ("Plots for clone_prev_MAE")
model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
ourcolors <- c("red", "blue", "green", "purple", "cyan")
label <- c("Epiclomal","EuclideanClust","DensityCut","HammingClust","PearsonClust")
list[big_df, crash] <- collect_data(model, label, ourcolors, criterion, "clone_prev_MAE")
plot_data(big_df, crash, "clone_prev_MAE")

##################
### plots V-measure ####
##################
### V-measure
print ("Plots for V-measure")
model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
ourcolors <- c("red", "blue", "green", "purple", "cyan")
label <- c("Epiclomal","EuclideanClust","DensityCut","HammingClust","PearsonClust")
list[big_df, crash] <- collect_data (model, label, ourcolors, criterion, "Vmeasure")
plot_data(big_df, crash, "Vmeasure")

##################
### plots nclusters ####
##################
print ("Plots for nclusters")
model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
ourcolors <- c("red", "blue", "green", "purple", "cyan")
label <- c("Epiclomal","EuclideanClust","DensityCut","HammingClust","PearsonClust")
list[big_df, crash] <- collect_data (model, label, ourcolors, criterion, "nclusters")
plot_data(big_df, crash, "nclusters")
