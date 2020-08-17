
# plotting Vmeasure, Hamming Distance, number of clusters, cluster prevalence mean absolute error and mean squared error

suppressMessages(library("argparse"))

library(stringr)


# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

#parser$add_argument("--summary_file", type="character", help="File path to the summary table") 
parser$add_argument("--input_dir", type="character", help="Input dir with the synthetic runs") 

parser$add_argument("--output_dir", type="character", default="output", help="Entire or just the beginning of the output directory file ")

parser$add_argument("--var", type="character", help="The variable value")

parser$add_argument("--criterion", type="character", help="The selection criterion for the best run: lower_bound or log_posterior")

args <- parser$parse_args() 

print(args)

#summary_table_file <- args$summary_file
#summary_table <- read.table(summary_table_file,header=TRUE,na.strings="NA",sep="\t")
#print(summary_table)
#print(class(summary_table))
#print(str(summary_table))

inputdir <- args$input_dir
outdir <- args$output_dir

criterion <- args$criterion

dir.create(outdir, showWarnings = FALSE)

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


collect_data <- function(model, number_data_sets, inputdir, criterion, measure_name){
# measure_name can be HD, Vmeasure, nclusters, cp_error

    #variable <- as.character(summary_table[,1])
    # variable is the value of the changed variable, for example if we are varying misspb, variable is 0.5, 0.6, 0.7, 0.8, 0.9, 0.95
    
    # MA 24 July 2020, now getting variable from the input dir directly
    subdirs <- list.files(inputdir, pattern = "*", recursive = FALSE)
    variable <- NULL
    for (sdir in subdirs) {
      variable <- c(variable, unlist(strsplit(sdir,"_"))[1])
    }
    variable <- unique(variable)
    

    if (measure_name == "HD") {
        measure_title <- "hamming distance"
        column <- "mean"
        fname <- "hdist"
    } else if (measure_name == "Vmeasure") {
        measure_title <- " V-measure"
        column <- "best_vmeasure"
        fname <- "results"
    } else if (measure_name == "nclusters") {
        measure_title <- " number of predicted clusters"
        column <- "nclusters_pred"
        fname <- "results"
    } else if (measure_name == "clone_prev_MAE") {
        measure_title <- " clone prevalence mean absolute error"
        column <- "clone_prev_MAE"
        fname <- "results"
    } else if (measure_name == "clone_prev_MSE") {
        measure_title <- " clone prevalence mean squared error"
        column <- "clone_prev_MSE"
        fname <- "results"
    } else if (measure_name == "uncertainty") {
        measure_title <- " uncertainty true positive rate"
        column <- "uncertainty"
        fname <- "results"
    }      
     
    if(grepl("MISSPB", args$var, fixed=TRUE))   {
        xlabel <- "Missing data proportion"
    } else if (grepl("NREGIONS", args$var, fixed=TRUE))   {
        xlabel <- "Number of regions"
        # xlabel <- "Percentage of loci different between clusters"
        # variable <- 100.0/variable
    } else if (grepl("CLONE_PREV", args$var, fixed=TRUE))   {
        xlabel <- "Clone prevalence"    
    } else if (grepl("ERROR", args$var, fixed=TRUE))   {
        xlabel <- "Error"
    } else if (grepl("NCELLS", args$var, fixed=TRUE))   {
        xlabel <- "Number of cells"    
    } else if (grepl("NCLONES", args$var, fixed=TRUE))   {
        xlabel <- "Number of clones"    
    } else if (grepl("NLOCI", args$var, fixed=TRUE))   {
        xlabel <- "Number of loci"   
    } else if (grepl("READSIZE", args$var, fixed=TRUE))   {
        xlabel <- "Read size"     
    } else if (grepl("CELL_CELL_VAR", args$var, fixed=TRUE))   {    
        xlabel <- "Cell to cell variability"
    } else if (grepl("PROP_CPG_FLIP", args$var, fixed=TRUE))   {
        xlabel <- "Proportions of CpGs flipped in the different region"     
    } else if (grepl("ndiffreg", args$var, fixed=TRUE))   {
      xlabel <- "Number of different regions"     
    } else if (grepl("subsamp_missp", args$var, fixed=TRUE))   {
      xlabel <- "Missing data proportion (subsampled)"  
    } else if (grepl("subsamp_ncells", args$var, fixed=TRUE))   {
      xlabel <- "Proportion of cells (subsampled)"  
    }
 
    
     
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
                print(paste0("Model ", model[m], " value ", variable[j], " number of data sets ", number_data_sets))
                for(i in 1:number_data_sets){
                    if (model[m] == "HammingClust" || model[m] == "DensityCut" || model[m] == "PearsonClust" || model[m] == "EuclideanClust") {
                        results_file <- paste0(inputdir, "/", variable[j], "_", i, "/simple_hclust/results_", model[m], ".txt")
                    } else if (model[m] == "region_bulk") {     
                        results_file <- paste0(inputdir, "/", variable[j], "_", i, "/result_region/",criterion,"/all_",fname,"_bestrun_region.tsv")
                    } else {     
                        results_file <- paste0(inputdir, "/", variable[j], "_", i, "/result_",model[m],"/",criterion,"/all_",fname,"_bestrun_",model[m],".tsv")
                    }                      
                    print (paste0('file is ', results_file))
                    t <- try(read.table(file=results_file,sep="\t",header=TRUE))   
                    if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
                        print("can't find file")
                        crash <- c(crash,0)
                        measure <- c(measure,NA)
                        VAR <- c(VAR,variable[j])
                        method <- c(method,model[m])                         
                    } else {
                        # for HD, we have to check 3 files
                        for (index in c(1,2,3)) {
							              t <- try(read.table(file=results_file,sep="\t",header=TRUE))   
							              if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
								              print(paste0("Can't find file ", results_file, " skipping it"))
								              next
							              }	
                            f <- read.table(file=results_file,sep="\t",header=TRUE)							
                            if (model[m] == "region_bulk") {
                                if (measure_name == "Vmeasure") {
                                    column <- "slsbulk_vmeasure"
                                }   
                                if (measure_name == "clone_prev_MAE") {
                                    column <- "slsbulk_clone_prev_MAE"
                                } 
                            } else {
                                if (measure_name == "Vmeasure") {
                                    column <- "best_vmeasure"
                                }   
                                if (measure_name == "clone_prev_MAE") {
                                    column <- "clone_prev_MAE"
                                }                                                
                            } 
                            new_method <- model[m]  
                            if (measure_name == "HD") {
                                if (index == 1) {
                                    original_file <- results_file
                                    new_method <- paste0 (new_method, "-uncorrected")
                                    results_file <- paste0(results_file,".corr.tsv")
                                } else if (index == 2) {
                                    new_method <- paste0 (new_method, "-corrected")
                                    results_file <- paste0(original_file,".naive.tsv")
                                } else if (index == 3) {
                                    new_method <- paste0 (new_method, "-naive")
                                }
                            } 
                            if (is.na(f[,column])) {       # this can happen for uncertainty if they, added on 9 Apr 2019
                                crash <- c(crash,0)
                            } else {
                                crash <- c(crash,1)
                            }    
                            measure <- c(measure,f[,column])
                            VAR <- c(VAR,variable[j])
                            method <- c(method,new_method) 
                            if (measure_name != "HD") {
                                break
                            }
                        }    
                    }      
                }
            }
        }

        VAR <- gsub("^1_", "", VAR)
        variable <- gsub("^1_", "", variable)

        big_df <- cbind(as.data.frame(measure),as.data.frame(crash),VAR,method)
        colnames(big_df) <- c("Measure","crash","VAR","method")
        str(big_df)
#label <- c("EpiclomalRegion","EpiclomalBulk","EpiclomalBasic","EuclideanClust","DensityCut","HammingClust","PearsonClust")
        #for (i in 1:length(model)){      
        #  big_df$method <- sub(pattern=model[i],x=big_df$method,replacement=label[i])      
        #}    
        #big_df$method <- factor(big_df$method,levels=label)

        # big_df$method <- factor(big_df$method,levels=model)
        big_df$VAR <- factor(big_df$VAR,levels=variable)

        # Now saving the data frame
        save(big_df, crash, file=savedfile)
    }  # end make the data files  

    print("Big DF")
    print(big_df)

    return(list("big_df"=big_df, "crash"=crash, "xlabel"=xlabel))
}    


##################
### box plots and line plots ####
##################
### hamming distance

#number_data_sets <- as.numeric(gsub("_|D", "", str_extract(summary_table_file, "_D[0-9]+_")))
number_data_sets <- 10

#initial_path_to_each_RUN <- paste0(unlist(strsplit(summary_table_file, "/FINAL"))[1],"/RUN/D_")
#print (initial_path_to_each_RUN)


all_regions_criterion <- "0_1_0.01"

#xlsfile <- paste0(outdir,"/SourceData.xlsx")
#if (file.exists(xlsfile)) {
#	file.remove(xlsfile)
#}

##################
### plots hamming distance ####
##################
print ("Plots for hamming distance")
#model <- c("region", "basic")
model <- c("region")
mylist <- collect_data (model, number_data_sets, inputdir, criterion, "HD")
#model <- c("region-uncorrected", "region-corrected", "region-naive", "basic-uncorrected", "basic-corrected", "basic-naive")
#model <- c("region-corrected", "region-naive", "basic-corrected", "basic-naive")
model <- c("region-uncorrected", "region-corrected", "region-naive")
plot_data(mylist$big_df, mylist$crash,  mylist$xlabel, model, "HD", add_points=FALSE)

#stop("Stopped")

##################
### plots clone_prev_MAE ####
##################
print ("Plots for clone_prev_MAE")
model <- c("region","region_bulk", "basic","EuclideanClust","DensityCut","HammingClust","PearsonClust")
# To add custom labels
# label <- c("EpiclomalRegion","EpiclomalBulk","EpiclomalBasic","EuclideanClust","DensityCut","HammingClust","PearsonClust")
mylist <- collect_data (model, number_data_sets, inputdir, criterion, "clone_prev_MAE")
plot_data(mylist$big_df, mylist$crash, mylist$xlabel, model, "clone_prev_MAE", add_points=FALSE)
#plot_data(mylist$big_df, mylist$crash, model, "clone_prev_MAE", add_points=FALSE)


##################
### plots V-measure ####
##################


### V-measure
print ("Plots for V-measure")
#model <- c("region","region_bulk", "basic","EuclideanClust","DensityCut","HammingClust","PearsonClust")
model <- c("region","basic","EuclideanClust","DensityCut","HammingClust","PearsonClust")
mylist <- collect_data (model, number_data_sets, inputdir, criterion, "Vmeasure")
plot_data(mylist$big_df, mylist$crash,  mylist$xlabel, model, "Vmeasure", add_points=FALSE)

##################
### plots nclusters ####
##################
print ("Plots for nclusters")
model <- c("region","basic","EuclideanClust","DensityCut","HammingClust","PearsonClust")
mylist <- collect_data (model, number_data_sets, inputdir, criterion, "nclusters")
plot_data(mylist$big_df, mylist$crash,  mylist$xlabel, model, "nclusters", add_points=FALSE)

##################
### plots uncertainty ####
##################

print ("Plots for uncertainty")
model <- c("region")
mylist <- collect_data (model, number_data_sets, inputdir, criterion, "uncertainty")
plot_data(mylist$big_df, mylist$crash,  mylist$xlabel, model, "uncertainty", add_points=FALSE)

