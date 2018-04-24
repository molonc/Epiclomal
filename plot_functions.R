
# plotting functions for the plot_final_results*.R files
##########################################

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##########################################

summarySE_new <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                          conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     median = median (xx[[col]], na.rm=na.rm),
                     first_quartile = as.vector(quantile (xx[[col]], probs=0.25, na.rm=na.rm)),
                     third_quartile = as.vector(quantile (xx[[col]], probs=0.75, na.rm=na.rm))
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# I tryed to write this function for the case when we have 10.00 clusters, but it doesn't work
reformat <- function(x) {
    rx <- 0
    if(max(x, na.rm=T)>=10) {
        rx <- format(round(x,1),nsmall=1)
    } else {
        rx <- format(round(x,2),nsmall=2)
    }
    return (rx)
}

##########################################

plot_data <- function(model, label, ourcolors, criterion, measure_name){
  # measure_name can be HD, Vmeasure, nclusters, cp_error
  
  variable <- datasets
  # variable is the value of the changed variable, for example if we are varying misspb, variable is 0.5, 0.6, 0.7, 0.8, 0.9, 0.95
  
  if (measure_name == "HD") {
    measure_title <- "Hamming Distance"
    column <- "mean"
    fname <- "hdist"
  } else if (measure_name == "Vmeasure") {
    measure_title <- "V-measure"
    column <- "best_vmeasure"
    fname <- "results"
  } else if (measure_name == "nclusters") {
    measure_title <- "Number of predicted clusters"
    column <- "nclusters_pred"
    fname <- "results"
  } else if (measure_name == "clone_prev_MAE") {
    measure_title <- "Clone prevalence MAE"
    column <- "clone_prev_MAE"
    fname <- "results"
  } else if (measure_name == "clone_prev_MSE") {
    measure_title <- "Clone prevalence MSE"
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
  
  
  
  # sub_big_df <- ddply(big_df, .(VAR,method),summarise,crash_perc=100*(1-mean(crash)))
  sub_big_df <- ddply(big_df, .(VAR,method),summarise,crash_perc=(1-mean(crash)))
  
  print(sub_big_df)
  
  # plot the box plots
  pHD <- ggplot(big_df, aes(x=method, y=Measure, fill=method)) +
    # The next line writes the y labels in the format x.xx so it aligns well with the bar plot.
    scale_y_continuous(labels = function(x) format(round(x,2),nsmall=2)) +
    geom_boxplot() + 
    #geom_boxplot(show.legend=F) + 
    facet_grid(~VAR) +
    #ggtitle(xlabel) +
    labs(x="", y = measure_title) 
  pHD <- pHD + 
    #guides(fill=FALSE) +
    theme(plot.title = element_text(size=20), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=20, colour= "black"), 
          # axis.text.x  = element_blank()
          axis.text.y  = element_text(size=16, colour= "black"),
          #panel.background = element_rect(fill="white",colour = 'black'), 
          axis.title.y =element_text(size=20), 
          axis.title.x=element_text(size=20),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position="none",
          legend.text=element_text(size=16)
          #strip.text.x = element_text(size =16)
    )
  pHD <- pHD + scale_fill_manual(values=ourcolors)  
  
  # ggsave(pHD,file=paste0(outdir,"/boxplot_",measure_name,"_",criterion,".pdf"),width=13.1,height=10.6)  
  
  
  # plot bar plots for crash
  bHD <-ggplot(sub_big_df, aes(x=method, y=crash_perc, fill=method)) +
    scale_y_continuous(breaks=seq(0,1,0.50), labels = function(x) format(round(x,2),nsmall=2)) +
    geom_bar(stat="identity") + facet_grid(~VAR) +
    # ggtitle(xlabel) +
    labs(x="", y = paste0("Unsuccess")) 
  bHD <- bHD + theme(plot.title = element_text(size=20), 
                     #axis.text.x  = element_text(angle=90, vjust=0.5, size=16, colour= "black"), 
                     axis.text.x  = element_blank(),
                     axis.text.y  = element_text(size=16, colour= "black"),
                     #panel.background = element_rect(fill="white",colour = 'black'), 
                     axis.title.y =element_text(size=20), 
                     axis.title.x=element_text(size=20),
                     legend.position="top",
                     legend.text=element_text(size=16),
                     strip.text.x = element_text(size =20) )
                     
  bHD <- bHD +  scale_fill_manual(values=ourcolors)                   
  
  # ggsave(bHD,file=paste0(outdir,"/barplot_",measure_name,"_",criterion,".pdf"),width=13.1,height=10.6)  
  
  pdf(file=paste0(outdir,"/boxplot_",measure_name,"_",criterion,".pdf"),onefile=TRUE,width=13.1,height=10.6)
  grid.arrange(arrangeGrob(bHD,nrow=1,ncol=1), arrangeGrob(pHD,nrow=1,ncol=1),heights=c(2.9,10.6))
  dev.off()
  
  # print("using ggarrange")
  #
  # figure <- ggarrange(bHD, pHD,labels=NULL,
  #                     ncol = 1, nrow = 2) ## this function ggarrange may be useful one day
  #
  # print("saving plot done with ggarrange")
  #
  # ggsave(figure,file=paste0(outdir,"/TEST.pdf"),width=13.1,height=10.6) 
  #
  # stop()
  
  # plot the mean and median line plots
  #aggre <- c("mean")
  # TODO For some reason, it doesn't work for median, it says "need numeric data"
  aggre <- c("mean", "median")
  
  print(str(big_df))  ### Measure is numeric!
  #big_df$Measure <- as.numeric(as.character(big_df$Measure))
  
  agg_df <- summarySE_new(big_df, measurevar="Measure", groupvars=c("VAR","method"),na.rm=TRUE)
  print(agg_df)
  
  for (agg in aggre) {
    
    # agg_df <- aggregate(big_df, by=list(Method=big_df$method, VAR=big_df$VAR), FUN=agg, na.rm=FALSE)   
    # print(paste0(agg, " DF"))
    # print(agg_df)    
    
    # agg_df <- summarySE(big_df, measurevar="Measure", groupvars=c("VAR","method"),na.rm=TRUE)
    # print(agg_df)
    
    # agg_df <- summarySE_new(big_df, measurevar="Measure", groupvars=c("VAR","method"),na.rm=TRUE)
    # print(agg_df)
    # 
    print(agg)
    
    
    
    # The errorbars overlapped, so use position_dodge to move them horizontally
    pd <- position_dodge(0.1) # move them .05 to the left and right, if 0 no move happens
    
    if (agg == "mean" ){
      
      pHD <- ggplot(agg_df, aes(x=VAR, y=Measure, group=method,colour=method)) +
        geom_errorbar(aes(ymin=Measure-se, ymax=Measure+se), width=.1,position=pd) +
        geom_line(aes(color=method), size=3,position=pd) + 
        geom_point(position=pd,size=4) +
        labs(x=xlabel, y = paste0(measure_title, " (", agg, ")")) 
      pHD <- pHD + theme(axis.text.y  = element_text(size=20, colour= "black"),
                         axis.text.x  = element_text(size=20, colour= "black"),      
                         axis.title.y =element_text(size=22), 
                         axis.title.x=element_text(size=22),
                         legend.text=element_text(size=16),
                         legend.position="top",
                         strip.text.x = element_text(size =16) )  
      pHD <- pHD + scale_color_manual(values=ourcolors)                           
      
      
    }
    
    if(agg == "median"){
      
      # The black error bars
      # pHD <- ggplot(agg_df, aes(x=VAR, y=median, group=method)) +
      #  geom_errorbar(aes(ymin=(median-(median-first_quartile)), ymax=(median+(third_quartile-median))), width=.1) +
      #   geom_line(aes(color=method), size=3) + 
      #   #geom_point() +
      #   labs(x=xlabel, y = paste0(agg, " ", measure_title))       
      
      pHD <- ggplot(agg_df, aes(x=VAR, y=median, group=method,colour=method)) +
        geom_errorbar(aes(ymin=Measure-se, ymax=Measure+se), width=.1,position=pd) +
        geom_line(aes(color=method), size=3,position=pd) + 
        geom_point(position=pd,size=4) +
        labs(x=xlabel, y = paste0(measure_title, " (", agg, ")")) 
      pHD <- pHD + theme(axis.text.y  = element_text(size=20, colour= "black"),
                         axis.text.x  = element_text(size=20, colour= "black"),
                         axis.title.y =element_text(size=22), 
                         axis.title.x=element_text(size=22),
                         legend.text=element_text(size=16),
                         legend.position="top",
                         strip.text.x = element_text(size =16) )    
      pHD <- pHD + scale_color_manual(values=ourcolors)
      
    }
    
    ggsave(pHD,file=paste0(outdir,"/lineplot_", agg, "_",measure_name,"_",criterion,".pdf"),width=13,height=10)              
  }
  
}