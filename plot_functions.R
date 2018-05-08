
library(ggplot2)
library(gridExtra)
library(plyr)


all_regions_criterion <- "0_1_0.01"


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

set_colors_and_labels <- function(model) {
  ourcolors <- vector(mode = "character", length = length(model))
  label <- vector(mode = "character", length = length(model))
  for (m in 1:length(model)) {
    if (model[m] == "region") {
      ourcolors[m] <- "red"
      label[m] <- "EpiclomalRegion"
    } else if (model[m] == "region_bulk") {
      ourcolors[m] <- "orange"
      label[m] <- "EpiclomalBulk"            
    } else if (model[m] == "basic") {
      ourcolors[m] <- "yellow"
      label[m] <- "EpiclomalBasic"            
    } else if (model[m] == "EuclideanClust" || model[m] == "Hclust") {
      ourcolors[m] <- "green"
      label[m] <- "EuclideanClust"            
    } else if (model[m] == "DensityCut" || model[m] == "densitycut") {
      ourcolors[m] <- "skyblue"
      label[m] <- "DensityCut"               
    } else if (model[m] == "HammingClust" || model[m] == "PBALclust") {
      ourcolors[m] <- "royalblue3"
      label[m] <- "HammingClust"               
    } else if (model[m] == "PearsonClust" || model[m] == "Pearsonclust") {
      ourcolors[m] <- "purple"
      label[m] <- "PearsonClust"               
    }        
  }
  print(model)
  print(ourcolors)
  print(label)
  return(list("colors"=ourcolors,"labels"=label))
}

########################################## 


multi_match <- function(x, table){
  # returns initial indicies of all substrings in table which match x
  if(length(table) < length(x)){
    return(NA)
  }else{
    check_mat <- matrix(nrow = length(x), ncol = length(table))
    for(i in 1:length(x)){
      check_mat[i,] <- table %in% x[i]
    }
    out <- vector(length = length(table))
    for(i in 1:(length(table)-(length(x)-1))){
      check <- vector(length=length(x))
      for(j in 1:length(x)){
        check[j] <- check_mat[j,(i+(j-1))]
      }
      out[i] <- all(check)
    }
    if(length(which(out))==0){
      return(NA)
    }else{
      return(which(out))
    }
  }
}


##############################################

grab_point <- function(x,big_df,criterion){
  
  #x_tmp <- big_df[multi_match(x,big_df$Measure):(multi_match(x,big_df$Measure)+length(x)-1),]
  
  #if NOT including that particular criterion in the box plot do the following then:
  
  #print(length(x))
  
  x_tmp <- big_df[multi_match(x,big_df$Measure):(multi_match(x,big_df$Measure)+length(x)),]
  
  print("camila")
  print(x_tmp)
  print(dim(x_tmp))
  
  if(sum(x_tmp$replicate == criterion) != 0){
    y <- x_tmp[which(x_tmp$replicate == criterion),]$Measure 
  }else{
    y <- NA    
  }
  
  print(y)
  return(y)  
}

################################################

plot_data <- function(big_df, crash, model, measure_name) {  
  
  our <- set_colors_and_labels(model)    
  ourcolors <- our$colors
  label <- our$label
  
  ## changing variable names    
  for (i in 1:length(model)){
    big_df$method <- sub(pattern=paste0("^",model[i],"$"),x=big_df$method,replacement=label[i])
  }    
  big_df$method <- factor(big_df$method,levels=label)
  #big_df$method <- factor(big_df$method,levels=method)
  
  if (measure_name == "HD") {
    measure_title <- "Hamming Distance"
    fname <- "hdist"
  } else if (measure_name == "Vmeasure") {
    measure_title <- "V-measure"
    fname <- "results"
  } else if (measure_name == "nclusters") {
    measure_title <- "Number of predicted clusters"
    fname <- "results"
  } else if (measure_name == "clone_prev_MAE") {
    measure_title <- "Clone prevalence MAE"
    fname <- "results"
  } else if (measure_name == "clone_prev_MSE") {
    measure_title <- "Clone prevalence MSE"
    fname <- "results"
  } else if (measure_name == "uncertainty") {
    measure_title <- "Uncertainty true positive rate"
    fname <- "results"
  }      
  
  xlabel = "Data set"  
  # sub_big_df <- ddply(big_df, .(VAR,method),summarise,crash_perc=100*(1-mean(crash)))
  sub_big_df <- ddply(big_df, .(VAR,method),summarise,crash_perc=(1-mean(crash)))   # big_df[['crash']])))
  
  print(sub_big_df)
  
  big_df_s <- subset(big_df, replicate != all_regions_criterion)
  
  sub_big_df <- ddply(big_df_s, .(VAR,method),summarise,crash_perc=(1-mean(crash)))   # big_df[['crash']])))
  
  print(sub_big_df)
  
  # print("testing")
  # print(big_df)
  # print(str(big_df))
  # 
  # sub <- subset(big_df,VAR == "InHouse" & method == "EpiclomalRegion")
  # print(sub)
  # 
  # x <- sub$Measure
  # x
  # 
  # print(big_df$Measure %in% x)
  # 
  # print(multi_match(x,big_df$Measure))
  # 
  # print(big_df[multi_match(x,big_df$Measure):(multi_match(x,big_df$Measure)+length(x)-1),])
  # 
  # x_tmp <- big_df[multi_match(x,big_df$Measure):(multi_match(x,big_df$Measure)+length(x)-1),]
  # 
  # criterion <-  "0_0.95_10000"
  # 
  # print(x_tmp[which(x_tmp$replicate == criterion),]$Measure) 
  
  # big_df$Measure[24] <- 0.25
  
  
  
  ### plot the box plots
  pHD <- ggplot(big_df_s, aes(x=method, y=Measure, fill=method)) +
    #pHD <- ggplot(big_df, aes(x=method, y=Measure, fill=method)) +
    # The next line writes the y labels in the format x.xx so it aligns well with the bar plot.
    scale_y_continuous(labels = function(x) format(round(x,2),nsmall=2)) +
    geom_boxplot() + 
    stat_summary(fun.y=grab_point,fun.args=list(big_df,criterion=all_regions_criterion), geom = "point",position=position_dodge(width=0.75),size=4,colour="blue",shape=18) +
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
  
  pdf(file=paste0(outdir,"/boxplot_",measure_name,".pdf"),onefile=TRUE,width=13.1,height=10.6)
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
    
    ggsave(pHD,file=paste0(outdir,"/lineplot_", agg, "_",measure_name,".pdf"),width=13,height=10)              
  }
  
}