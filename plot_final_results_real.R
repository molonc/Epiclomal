
# plotting Vmeasure, Hamming Distance, number of clusters, cluster prevalence mean absolute error and mean squared error

suppressMessages(library("argparse"))
library(ggplot2)
library(gridExtra)
library(stringr)
library(plyr)
#library(ggpubr)

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

parser$add_argument("--output_dir", type="character", default="output", help="Entire or just the beginning of the output directory file ")

parser$add_argument("--criterion", type="character", help="The selection criterion for the best run: DIC_measure_gainthr0.05 DIC_measure_gainthr0.1")

datasets <- c("Aparicio",
    "Smallwood2014",
    "Hou2016",
    "Luo2017",
    "Farlik2016")
datapaths <- c("../EPI-112_inhouse_data/FINAL_RESULTS",
    "../EPI-70_Smallwood2014/FINAL_RESULTS",
    "../EPI-105_scTrio/FINAL_RESULTS",
    "../EPI-106_Luo2017/FINAL_RESULTS",
    "../EPI-89_Farlik2016_all_union/FINAL_RESULTS")
simplepaths <- c("../EPI-112_inhouse_data/OUTPUT_epiclomal_INHOUSE/RUN/epiclomal_INHOUSE_",
    "../EPI-70_Smallwood2014/OUTPUT_epiclomal_Smallwood2014/RUN/epiclomal_Smallwood2014_",
    "../EPI-105_scTrio/OUTPUT_epiclomal_scTrio/RUN/epiclomal_scTrio_",
    "../EPI-106_Luo2017/OUTPUT_epiclomal_Luo2017_genebodies_500_clean_random_cells/RUN/epiclomal_Luo2017_genebodies_500_clean_random_cells_",
    "../EPI-89_Farlik2016_all_union/OUTPUT_epiclomal_Farlik2016_all_union/RUN/epiclomal_Farlik2016_all_union_")


# each replicate file should be in inputs, for example inputs/Smallwood2014_replicates.txt

args <- parser$parse_args() 

print(args)

outdir <- args$output_dir
dir.create(outdir, showWarnings = FALSE)

criterion <- args$criterion
dir.create(paste0(outdir,"/",criterion), showWarnings = FALSE)
outdir <- paste0(outdir,"/",criterion)

##################
### box plots ####
##################

plot_data <- function(model,label,criterion, measure_name){
# measure_name can be HD, Vmeasure, nclusters, cp_error

    variable <- datasets
    # variable is the value of the changed variable, for example if we are varying misspb, variable is 0.5, 0.6, 0.7, 0.8, 0.9, 0.95

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
        save(big_df, file=savedfile)
 }  # end make the data files  
  
  
    
  sub_big_df <- ddply(big_df, .(VAR,method),summarise,crash_perc=100*(1-mean(crash)))

  print(sub_big_df)
  
    # plot the box plots
    pHD <- ggplot(big_df, aes(x=method, y=Measure, fill=method)) +
      geom_boxplot() + 
      #geom_boxplot(show.legend=F) + 
      facet_grid(~VAR) +
      #ggtitle(xlabel) +
      labs(x="", y = paste0("Cell-based ", measure_title)) 
      pHD <- pHD + 
            #guides(fill=FALSE) +
            theme(plot.title = element_text(size=20), 
            axis.text.x  = element_text(angle=90, vjust=0.5, size=16, colour= "black"), 
            # axis.text.x  = element_blank()
            axis.text.y  = element_text(size=16, colour= "black"),
            #panel.background = element_rect(fill="white",colour = 'black'), 
            axis.title.y =element_text(size=16), 
            axis.title.x=element_text(size=20),
            strip.background = element_blank(),
            strip.text.x = element_blank()
            #legend.position="none",
            #strip.text.x = element_text(size =16)
            )

# ggsave(pHD,file=paste0(outdir,"/boxplot_",measure_name,"_",criterion,".pdf"),width=13.1,height=10.6)  
     

    # plot bar plots for crash
    bHD <-ggplot(sub_big_df, aes(x=method, y=crash_perc, fill=method)) +
      geom_bar(stat="identity") + facet_grid(~VAR) +
      ggtitle(xlabel) +
      labs(x="", y = paste0("Unsuccessful runs %")) 
    bHD <- bHD + theme(plot.title = element_text(size=20), 
                   #axis.text.x  = element_text(angle=90, vjust=0.5, size=16, colour= "black"), 
                   axis.text.x  = element_blank(),
                   axis.text.y  = element_text(size=10, colour= "black"),
                   #panel.background = element_rect(fill="white",colour = 'black'), 
                   axis.title.y =element_text(size=12), 
                   axis.title.x=element_text(size=20),
                   strip.text.x = element_text(size =16) )

# ggsave(bHD,file=paste0(outdir,"/barplot_",measure_name,"_",criterion,".pdf"),width=13.1,height=10.6)  

pdf(file=paste0(outdir,"/boxplot_",measure_name,"_",criterion,".pdf"),onefile=TRUE,width=13.1,height=10.6)
grid.arrange(arrangeGrob(bHD,nrow=1,ncol=1), arrangeGrob(pHD,nrow=1,ncol=1),heights=c(2.5,10.6))
dev.off()

#figure <- ggarrange(bHD, pHD,
#                    ncol = 1, nrow = 2) ## this function ggarrange may be useful one day

    
    # plot the mean and median line plots
    aggre <- c("mean")
    # TODO For some reason, it doesn't work for median, it says "need numeric data"
    # aggre <- c("mean", "median")
    big_df$Measure <- as.numeric(as.character(big_df$Measure))
    for (agg in aggre) {
        agg_df <- aggregate(big_df, by=list(Method=big_df$method, VAR=big_df$VAR), FUN=agg, na.rm=FALSE)    
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
model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
label <- c("Epiclomal","EuclideanClust","DensityCut","HammingClust","PearsonClust")
plot_data (model, label,criterion, "clone_prev_MAE")



##################
### plots V-measure ####
##################
### V-measure
print ("Plots for V-measure")
model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
label <- c("Epiclomal","EuclideanClust","DensityCut","HammingClust","PearsonClust")
plot_data (model,label, criterion, "Vmeasure")


##################
### plots nclusters ####
##################
print ("Plots for nclusters")
model <- c("region", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
label <- c("Epiclomal","EuclideanClust","DensityCut","HammingClust","PearsonClust")
plot_data (model, label,criterion, "nclusters")
