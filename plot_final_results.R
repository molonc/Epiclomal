
# plotting Vmeasure, Hamming Distance, number of clusters, cluster prevalence mean absolute error and mean squared error

suppressMessages(library("argparse"))
library(ggplot2)
library(stringr)

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

parser$add_argument("--summary_file", type="character", help="File path to the summary table") 

parser$add_argument("--output_dir", type="character", default="output", help="Entire or just the beginning of the output directory file ")

parser$add_argument("--var", type="character", help="The variable value")

parser$add_argument("--criterion", type="character", help="The selection criterion for the best run: lower_bound or log_posterior")

args <- parser$parse_args() 

print(args)

summary_table_file <- args$summary_file

summary_table <- read.table(summary_table_file,header=TRUE,na.strings="NA",sep="\t")

print(summary_table)

print(class(summary_table))

print(str(summary_table))

outdir <- args$output_dir

criterion <- args$criterion


##################
### box plots ####
##################

plot_data <- function(model, number_data_sets, initial_path_to_each_RUN, summary_table, criterion, measure_name){
# measure_name can be HD, Vmeasure, nclusters, cp_error

    variable <- as.character(summary_table[,1])
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
     
    if(grepl("MISSPB", args$var, fixed=TRUE))   {
        xlabel <- "Missing data proportion"
    } else if (grepl("NREGIONS", args$var, fixed=TRUE))   {
        xlabel <- "Percentage of loci different between clusters"
        x <- 100.0/x
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

        counts <- c(0,0)

        for(m in 1:length(model)){
            for(j in 1:length(variable)){
                print(paste0("Model ", model[m], " value ", variable[j], " number of data sets ", number_data_sets))
                for(i in 1:number_data_sets){
                    if (model[m] == "PBALclust" || model[m] == "densitycut" || model[m] == "Pearsonclust" || model[m] == "Hclust") {
                        results_file <- paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/simple_hclust/results_", model[m], ".txt")
                    } else {     
                        results_file <- paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/results_",model[m],"/",criterion,"/all_",fname,"_bestrun_",model[m],".tsv")
                    }    
                    # print (paste0('file is ', results_file))
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

        #print("Big DF")
        #print(big_df)
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
    # big_df$Measure <- as.numeric(as.character(big_df$Measure))
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
### hamming distance

number_data_sets <- as.numeric(gsub("_|D", "", str_extract(summary_table_file, "_D[0-9]+_")))


initial_path_to_each_RUN <- paste0(unlist(strsplit(summary_table_file, "/FINAL"))[1],"/RUN/D_")
print (initial_path_to_each_RUN)


##################
### plots clone_prev_MAE ####
##################
print ("Plots for clone_prev_MAE")
# model <- c("basic","basic_munok","region","region_munok","PBALclust","densitycut")
model <- c("region", "basic", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
plot_data (model, number_data_sets, initial_path_to_each_RUN, summary_table, criterion, "clone_prev_MAE")

##################
### plots clone_prev_MSE ####
##################
# print ("Boxplots for clone_prev_MSE")
# model <- c("basic","basic_munok","region","region_munok","PBALclust","densitycut")
# plot_data (model, number_data_sets, initial_path_to_each_RUN, summary_table, criterion, "clone_prev_MSE")

##################
### plots hamming distance ####
##################
print ("Plots for hamming distance")
# model <- c("basic","basic_munok","region","region_munok")
model <- c("region", "basic")
plot_data (model, number_data_sets, initial_path_to_each_RUN, summary_table, criterion, "HD")

##################
### plots V-measure ####
##################
### V-measure
print ("Plots for V-measure")
# 
# model <- c("basic","basic_munok","region","region_munok","PBALclust","densitycut")
model <- c("region", "basic", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
plot_data (model, number_data_sets, initial_path_to_each_RUN, summary_table, criterion, "Vmeasure")

##################
### box plots V-measure ####
##################
### V-measure
#print ("Boxplots for V-measure")
#model <- c("basic","region")
#boxplot_data (model, number_data_sets, initial_path_to_each_RUN, summary_table, criterion, "Vmeasure")

##################
### plots nclusters ####
##################
print ("Plots for nclusters")
# model <- c("basic","basic_munok","region","region_munok","PBALclust","densitycut")
model <- c("region", "basic", "Hclust", "densitycut", "PBALclust", "Pearsonclust")
plot_data (model, number_data_sets, initial_path_to_each_RUN, summary_table, criterion, "nclusters")

# ###############################
# #### Code for plotting ELBo ###
# ###############################
# 
# path_to_elbo <- "/Users/cdesouza/Documents/shahlab15/mandronescu/EPI-73_run_synthetic_pipeline/OUTPUT_MISSPB_10CLONES/RUN/D_MISSPB_0.9_1_epiclomal_synthetic"
# 
# 
# big_df <- NULL
# 
# initials <- seq(0:49) -1
# 
# method <- NULL
# 
# for( i in 1:length(initials)){
#   
#   print(i)
#   print(initials[i])
#   
#   if(initials[i] >= 0 & initials[i] <= 9 ){
# 
#   hclust_crash <- as.numeric(read.table(file=paste0(path_to_elbo,"/outputs/simple_hclust/hclust_region_crash.tsv")))
# 
#   if(hclust_crash == 0){
#     
#   file <- paste0(path_to_elbo,"/outputs/epiclomal_basic/",initials[i],"/params.yaml")
#   
#   lines <- readLines(file)
# 
#   converged <- sum(sapply(sapply(lines, grep, pattern = "converged: true", value = TRUE), length) == 1)
#   
#   if(converged == 1){
#   
#   tmp <- as.numeric(sub("lower_bound: ", "", sapply(lines, grep, pattern = "lower_bound", value = TRUE)))
#   elbo <- tmp[!is.na(tmp)]
#   rm(tmp)
#   
#   tmp <- as.numeric(sub("CPU_time_seconds: ", "", sapply(lines, grep, pattern = "CPU_time_seconds", value = TRUE)))
#   cpu_time <- tmp[!is.na(tmp)]
#   rm(tmp)
#   
#   tmp <- as.numeric(sub("Vmeasure: ", "", sapply(lines, grep, pattern = "Vmeasure", value = TRUE)))
#   vmeasure <- tmp[!is.na(tmp)]
#   rm(tmp)
#   
#   big_df <- rbind(big_df,c(elbo,cpu_time,vmeasure))
#   
#   method <- c(method,"hclust")
#   
#   } 
#  
# }
# 
#   }
#   
#   if(initials[i] >= 10 & initials[i] <= 19 ){
#     
#     pbal_crash <- as.numeric(read.table(file=paste0(path_to_elbo,"/outputs/simple_hclust/PBAL_crash.tsv")))
#     
#     if(pbal_crash == 0){
#       
#       file <- paste0(path_to_elbo,"/outputs/epiclomal_basic/",initials[i],"/params.yaml")
#       
#       lines <- readLines(file)
#       
#       converged <- sum(sapply(sapply(lines, grep, pattern = "converged: true", value = TRUE), length) == 1)
#       
#       if(converged == 1){
#         
#         tmp <- as.numeric(sub("lower_bound: ", "", sapply(lines, grep, pattern = "lower_bound", value = TRUE)))
#         elbo <- tmp[!is.na(tmp)]
#         rm(tmp)
#         
#         tmp <- as.numeric(sub("CPU_time_seconds: ", "", sapply(lines, grep, pattern = "CPU_time_seconds", value = TRUE)))
#         cpu_time <- tmp[!is.na(tmp)]
#         rm(tmp)
#         
#         tmp <- as.numeric(sub("Vmeasure: ", "", sapply(lines, grep, pattern = "Vmeasure", value = TRUE)))
#         vmeasure <- tmp[!is.na(tmp)]
#         rm(tmp)
#         
#         big_df <- rbind(big_df,c(elbo,cpu_time,vmeasure))
#         
#         method <- c(method,"pbal")
#         
#       } 
#       
#     }
#     
#   }
#   
#   
#   if(initials[i] > 19 ){
#     
#     
#       
#       file <- paste0(path_to_elbo,"/outputs/epiclomal_basic/",initials[i],"/params.yaml")
#       
#       lines <- readLines(file)
#       
#       converged <- sum(sapply(sapply(lines, grep, pattern = "converged: true", value = TRUE), length) == 1)
#       
#       if(converged == 1){
#         
#         tmp <- as.numeric(sub("lower_bound: ", "", sapply(lines, grep, pattern = "lower_bound", value = TRUE)))
#         elbo <- tmp[!is.na(tmp)]
#         rm(tmp)
#         
#         tmp <- as.numeric(sub("CPU_time_seconds: ", "", sapply(lines, grep, pattern = "CPU_time_seconds", value = TRUE)))
#         cpu_time <- tmp[!is.na(tmp)]
#         rm(tmp)
#         
#         tmp <- as.numeric(sub("Vmeasure: ", "", sapply(lines, grep, pattern = "Vmeasure", value = TRUE)))
#         vmeasure <- tmp[!is.na(tmp)]
#         rm(tmp)
#         
#         big_df <- rbind(big_df,c(elbo,cpu_time,vmeasure))
#         
#         method <- c(method,"random")
#         
#       } 
#       
#     }
# 
# }
# 
# colnames(big_df) <- c("elbo","cpu_time","vmeasure")
# big_df <- cbind(as.data.frame(big_df),method)
# 
# outdir <- "/Users/cdesouza/Documents/"
# 
# pelbo <- ggplot(big_df, aes(x=method, y=elbo,fill=method)) +
#   geom_boxplot() + 
#   labs(x="initialization method", y = "elbo") +
#   theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12, colour= "black"), axis.text.y  = element_text(size=15, colour= "black"),
#         #panel.background = element_rect(fill="white",colour = 'black'), 
#         axis.title.y =element_text(size=15), axis.title.x=element_text(size=15),
#         strip.text.x = element_text(size =12) )
# 
# ggsave(pelbo,file=paste0(outdir,"boxplot_elbo_initialization.pdf"),width=13.1,height=10.6)
# 
# pcpu <- ggplot(big_df, aes(x=method, y=cpu_time,fill=method)) +
#   geom_boxplot() + 
#   labs(x="initialization method", y = "cpu time") +
#   theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12, colour= "black"), axis.text.y  = element_text(size=15, colour= "black"),
#         #panel.background = element_rect(fill="white",colour = 'black'), 
#         axis.title.y =element_text(size=15), axis.title.x=element_text(size=15),
#         strip.text.x = element_text(size =12) )
# 
# ggsave(pcpu,file=paste0(outdir,"boxplot_cpu_initialization.pdf"),width=13.1,height=10.6)
# 
# 
# pv <- ggplot(big_df, aes(x=method, y=vmeasure,fill=method)) +
#   geom_boxplot() + 
#   labs(x="initialization method", y = "V measure") +
#   theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12, colour= "black"), axis.text.y  = element_text(size=15, colour= "black"),
#         #panel.background = element_rect(fill="white",colour = 'black'), 
#         axis.title.y =element_text(size=15), axis.title.x=element_text(size=15),
#         strip.text.x = element_text(size =12) )
# 
# ggsave(pv,file=paste0(outdir,"boxplot_vmeasure_initialization.pdf"),width=13.1,height=10.6)









