
# plotting Vmeasure, Hamming Distance, number of clusters, cluster prevalence mean absolute error and mean squared error

suppressMessages(library("argparse"))
library(ggplot2)
library(stringr)

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

parser$add_argument("--output_dir", type="character", default="output", help="Entire or just the beginning of the output directory file ")

parser$add_argument("--criterion", type="character", help="The selection criterion for the best run: DIC_measure_gainthr0.05 DIC_measure_gainthr0.1")

parser$add_argument("--datasets", type="character", help="The data sets, e.g. 0.8_0.85_0.9_0.95")

parser$add_argument("--nreplicates", type="character", help="The number of times we did subsampling for each case")

#datasets <- c("0.7", "0.75", "0.8", "0.85", "0.9", "0.95")

# now we have 1 path for all subsamples

parser$add_argument("--datapath", type="character", help="Path with the epiclomal results, e.g. FINAL_RESULTS/")
# datapath <- "FINAL_RESULTS/"

# each replicate file should be in inputs, for example inputs/Smallwood2014_replicates.txt

args <- parser$parse_args() 

print(args)

outdir <- args$output_dir
dir.create(outdir, showWarnings = FALSE)

criterion <- args$criterion
dir.create(paste0(outdir,"/",criterion), showWarnings = FALSE)
outdir <- paste0(outdir,"/",criterion)

datasets <- unlist(strsplit(args$datasets, split="_"))
datapath <- args$datapath
nreplicates <- args$nreplicates

##################
### box plots ####
##################

plot_data <- function(model, criterion, measure_name){
# measure_name can be HD, Vmeasure, nclusters, cp_error

    variable <- datasets
    # variable is the value of the changed variable, for example if we are varying misspb, variable is 0.5, 0.6, 0.7, 0.8, 0.9, 0.95


    if (measure_name == "HD") {
        measure_title <- "hamming distance"
        column <- "hd_mean"
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
    }       
     
    xlabel = "Data set"
         
    savedfile <- paste0(outdir,"/data_",measure_name,"_",criterion,".Rda")
    if (file.exists(savedfile)) {
        print("File already exists, loading it")
        load(savedfile)
    } else {
        print("File doesn't exist, creating it")
        method <- NULL
        VAR <- NULL
        measure <- NULL

        counts <- c(0,0)

        for(m in 1:length(model)){
            for(j in 1:length(variable)){
                print(paste0("Model ", model[m], " data set ", variable[j], " number of replicates ", nreplicates))
                for(i in 1:nreplicates){
                    if (model[m] == "HammingClust" || model[m] == "DensityCut" || model[m] == "PearsonClust" || model[m] == "EuclideanClust") {
                        results_file <- paste0(datapath,"/", variable[j], "_",i,"/simple_hclust/results_", model[m], ".txt")                        
                    } else {     
                        results_file <- paste0(datapath,"/",variable[j], "_", i,"/result_",model[m],"/",criterion,"/all_results_bestrun_",model[m],".tsv")
                    }    
                    print (paste0('Results file is ', results_file))
                    t <- try(read.table(file=results_file,sep="\t",header=TRUE))   
                    if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
                        print("can't find file")
                        counts[j] <- counts[j]+1 
                    } else {
                        f <- read.table(file=results_file,sep="\t",header=TRUE)
                        measure <- c(measure,f[,column])
                        VAR <- c(VAR,variable[j])
                        method <- c(method,model[m]) 
                    }      
                }
            }
        }


        big_df <- cbind(as.data.frame(measure),VAR,method)
        colnames(big_df) <- c("Measure","VAR","method")
        str(big_df)

        big_df$method <- factor(big_df$method,levels=model)
        big_df$VAR <- factor(big_df$VAR,levels=variable)

        print("Big DF")
        print(big_df)
        # Now saving the data frame
        save(big_df, file=savedfile)
    }  # end make the data files  

    # plot the box plots
    pHD <- ggplot(big_df, aes(x=method, y=Measure, fill=method)) +
      geom_boxplot() + facet_grid(~VAR) +
      ggtitle(xlabel) +
      labs(x="", y = paste0("Cell-based ", measure_title)) 
      pHD <- pHD + theme(plot.title = element_text(size=20), 
            axis.text.x  = element_text(angle=90, vjust=0.5, size=16, colour= "black"), 
            # axis.text.x  = element_blank()
            axis.text.y  = element_text(size=20, colour= "black"),
            #panel.background = element_rect(fill="white",colour = 'black'), 
            axis.title.y =element_text(size=20), 
            axis.title.x=element_text(size=20),
            strip.text.x = element_text(size =16) )
        

    ggsave(pHD,file=paste0(outdir,"/boxplot_",measure_name,"_",criterion,".pdf"),width=13.1,height=10.6)    
    
    # plot the mean and median line plots
    
    aggre <- c("mean")
    # TODO For some reason, it doesn't work for median, it says "need numeric data"
    # aggre <- c("mean", "median")
    big_df$Measure <- as.numeric(as.character(big_df$Measure))
    for (agg in aggre) {
        agg_df <- aggregate(big_df, by=list(Method=big_df$method, VAR=big_df$VAR), FUN=agg, na.rm=FALSE)[,1:3]    
        print(paste0(agg, " DF"))
        print(agg_df)    
        pHD <- ggplot(agg_df, aes(x=VAR, y=Measure, group=Method)) +
            geom_line(aes(color=Method), size=3) + 
            labs(x=xlabel, y = paste0(agg, " ", measure_title)) 
        pHD <- pHD + theme(axis.text.y  = element_text(size=20, colour= "black"),
            axis.title.y =element_text(size=20), axis.title.x=element_text(size=20),
            strip.text.x = element_text(size =16) )    

        ggsave(pHD,file=paste0(outdir,"/lineplot_", agg, "_",measure_name,"_",criterion,".pdf"),width=15,height=10)              
    }        
}



##################
### box plots and line plots ####
##################


##################
### plots clone_prev_MAE ####
##################
print ("Plots for clone_prev_MAE")
model <- c("region", "basic", "EuclideanClust", "DensityCut", "HammingClust", "PearsonClust")
plot_data (model, criterion, "clone_prev_MAE")

##################
### plots V-measure ####
##################
### V-measure
print ("Plots for V-measure")
model <- c("region", "basic", "EuclideanClust", "DensityCut", "HammingClust", "PearsonClust")
plot_data (model, criterion, "Vmeasure")

##################
### plots nclusters ####
##################
print ("Plots for nclusters")
model <- c("region", "basic", "EuclideanClust", "DensityCut", "HammingClust", "PearsonClust")
plot_data (model, criterion, "nclusters")


##################
### plots hamming distance ####
##################
#print ("Plots for HD")
#model <- c("region")
#plot_data (model, criterion, "HD")
