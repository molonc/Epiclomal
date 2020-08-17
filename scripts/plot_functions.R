
library(ggplot2)
library(gridExtra)
library(plyr)
#library(xlsx)


# one_filter = 1      # only one filter -- for the main text
 one_filter = 0      # all the filters -- for sup
# inputs directory has to correspond with the right filters




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
  # the colors from RcolorBrewer set1:
  # http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=8  
  ourcolors <- vector(mode = "character", length = length(model))
  label <- vector(mode = "character", length = length(model))
  for (m in 1:length(model)) {
    if (model[m] == "region") {
      ourcolors[m] <- "#e41a1c"        # was red, then red3
      label[m] <- "EpiclomalRegion"
    } else if (model[m] == "region-uncorrected") {
      ourcolors[m] <- "#a65628"
      label[m] <- "EpiclomalRegion-unadjusted"
    } else if (model[m] == "region-corrected") {
      ourcolors[m] <- "#e41a1c"
      label[m] <- "EpiclomalRegion-adjusted"      
    } else if (model[m] == "region-naive") {
      ourcolors[m] <- "grey15"     # brown
      label[m] <- "EpiclomalRegion-naive"      
    } else if (model[m] == "region_bulk") {
      ourcolors[m] <- "#ff7f00"      # was orange
      label[m] <- "EpiclomalBulk"            
    } else if (model[m] == "basic") {
      ourcolors[m] <- "gold"        # was yellow or gold
      label[m] <- "EpiclomalBasic"            
    } else if (model[m] == "basic-uncorrected") {
      ourcolors[m] <- "goldenrod"
      label[m] <- "EpiclomalBasic-unadjusted"
    } else if (model[m] == "basic-corrected") {
      ourcolors[m] <- "gold"
      label[m] <- "EpiclomalBasic-adjusted"     
    } else if (model[m] == "basic-naive") {
      ourcolors[m] <- "grey40"
      label[m] <- "EpiclomalBasic-naive"                  
    } else if (model[m] == "EuclideanClust" || model[m] == "Hclust") {
      ourcolors[m] <- "#4daf4a"       # was green or limegreen
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

########################################## 

grab_point <- function(x){
  
  df_tmp <- big_df[!is.na(big_df$Measure),]
  
  if(!is.na(multi_match(x,df_tmp$Measure)[1])){
    print("match found")
    
    x_tmp <- df_tmp[multi_match(x,df_tmp$Measure)[1]:(multi_match(x,df_tmp$Measure)[1]+length(x)-1),]
    #print(x_tmp)
    #print(as.character(x_tmp$VAR[1]))
    #print(as.character(x_tmp$method[1]))
    #print(big_df[which((big_df$VAR == as.character(x_tmp$VAR[1])) & (big_df$method == as.character(x_tmp$method[1]))),])
    #print(big_df[which((big_df$VAR == as.character(x_tmp$VAR[1])) & (big_df$method == as.character(x_tmp$method[1]))),1][1])
    
    return(big_df[which((big_df$VAR == as.character(x_tmp$VAR[1])) & (big_df$method == as.character(x_tmp$method[1]))),1][1])
    
  }else{
    print("no match")
    return(NA)
  }
  
}

################################################

plot_data <- function(big_df, crash, xlabel, model, measure_name, add_points) {  
  
  # select only what is in "model"
  big_df <- big_df[big_df$method %in% model,]
  print("Big df after selecting model")
  print(big_df)
  
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
  
  # xlabel = "Data set"  # Now xlabel is given
  
  #sub_big_df <- ddply(big_df, .(VAR,method),summarise,crash_perc=(1-mean(crash)))   # big_df[['crash']])))
  #print(sub_big_df)
  
  
  ### plot the box plots
  if(add_points==TRUE){
    
    big_df_s <- subset(big_df, replicate != all_regions_criterion)
    
    sub_big_df <- ddply(big_df_s, .(VAR,method),summarise,crash_perc=(1-mean(crash)))   # big_df[['crash']])))
    
    df_point <- big_df[which(big_df$replicate == all_regions_criterion),]
    
    df_point <- df_point[,c(1,3,4)]
    
    print(df_point)
    
    pHD <- ggplot(big_df_s, aes(x=method, y=Measure, fill=method)) +
      #pHD <- ggplot(big_df, aes(x=method, y=Measure, fill=method)) +
      # The next line writes the y labels in the format x.xx so it aligns well with the bar plot.
      scale_y_continuous(labels = function(x) format(round(x,2),nsmall=2)) +
      geom_boxplot() +    
      
      geom_point(data = df_point, aes(group = method), position = position_dodge(width = 0.75),
                 shape=18,color="blue",size=4) +
      
      #stat_summary(fun.y="mean", geom = "point",position=position_dodge(width=0.75),size=4,colour="blue",shape=18) +
      #stat_summary(fun.y=grab_point, geom = "point",position=position_dodge(width=0.75),size=4,colour="blue",shape=18) +
      
      #geom_boxplot(show.legend=F) + 
      facet_grid(~VAR) +
      labs(x="", y = measure_title)
  }
  
  if(add_points == FALSE){
    
    big_df_s <- big_df
    sub_big_df <- ddply(big_df_s, .(VAR,method),summarise,crash_perc=(1-mean(crash)))   # big_df[['crash']])))      
    pHD <- ggplot(big_df_s, aes(x=method, y=Measure, fill=method)) +
      #pHD <- ggplot(big_df, aes(x=method, y=Measure, fill=method)) +
      # The next line writes the y labels in the format x.xx so it aligns well with the bar plot.
      scale_y_continuous(labels = function(x) format(round(x,2),nsmall=2)) +
      geom_boxplot() +    
      facet_grid(~VAR) +
      labs(x="", y = measure_title)
    
  }  
  
  pHD <- pHD + 
    #guides(fill=FALSE) +
    theme(plot.title = element_text(size=25), 
          axis.text.x  = element_text(angle=90, vjust=0.5, size=25, colour= "black"), 
          # axis.text.x  = element_blank()
          axis.text.y  = element_text(size=21, colour= "black"),
          #panel.background = element_rect(fill="white",colour = 'black'), 
          axis.title.y =element_text(size=23), 
          axis.title.x=element_text(size=25),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position="none",
          legend.text=element_text(size=25)
          #strip.text.x = element_text(size =16)
    )
  pHD <- pHD + scale_fill_manual(values=ourcolors)  
  # plot bar plots for crash
  if (measure_name == "uncertainty") {
    ylabel <- "Undefined" 
  } else {
    ylabel <- "Failure"
  }    
  bHD <-ggplot(sub_big_df, aes(x=method, y=crash_perc, fill=method)) +
    scale_y_continuous(breaks=seq(0,1,0.50), limits=c(0,1), labels = function(x) format(round(x,2),nsmall=2)) +
    geom_bar(stat="identity") + facet_grid(~VAR) +
    # ggtitle(xlabel) +
    labs(x="", y = ylabel) 
  bHD <- bHD + theme(plot.title = element_text(size=25), 
                     #axis.text.x  = element_text(angle=90, vjust=0.5, size=16, colour= "black"), 
                     axis.text.x  = element_blank(),
                     axis.text.y  = element_text(size=21, colour= "black"),
                     #panel.background = element_rect(fill="white",colour = 'black'), 
                     axis.title.y =element_text(size=23), 
                     axis.title.x=element_text(size=25),
                     legend.position="top",
                     legend.text=element_text(size=22),
                     legend.title=element_text(size=25),
                     strip.text.x = element_text(size =23) )
  
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
  
  ######################################
  # plot the mean and median line plots
  #aggre <- c("mean")
  # TODO For some reason, it doesn't work for median, it says "need numeric data"
  aggre <- c("mean", "median")
  
  print(str(big_df))  ### Measure is numeric!
  #big_df$Measure <- as.numeric(as.character(big_df$Measure))
  
  agg_df <- summarySE_new(big_df, measurevar="Measure", groupvars=c("VAR","method"),na.rm=TRUE)
  
  #xlsfile <- paste0(outdir,"/SourceData.xlsx")
  #write.xlsx(agg_df, file=xlsfile, sheetName=measure_name,append=TRUE)
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
        
      if (measure_name == "HD") {
        lpos <- c(0.25,0.8)
        if (xlabel == "Number of loci")    {  lpos <- c(0.6,0.8)  }
        if (xlabel == "Clone prevalence")  {  lpos <- c(0.25,0.2) }
        if (xlabel == "Number of cells")   {  lpos <- c(0.45,0.8) }
        # lpos <- c(0.6,0.8)      # put legend on top left
		# for NCELLS lpos <- c(0.45,0.8)
		# for CLONE_PREV_500CELLS lpos <- c(0.25,0.2)
		# for NLOCI lpos <- c(0.6,0.8)
		# for others lpos <- c(0.25,0.8)		
      }  else {
        lpos <- "top"
      }
      print("Our colors")
      print(ourcolors)      
      pHD <- pHD + theme(axis.text.y  = element_text(size=25, colour= "black"),
                         axis.text.x  = element_text(size=25, colour= "black"),      
                         axis.title.y =element_text(size=27), 
                         axis.title.x=element_text(size=27),
                         legend.text=element_text(size=22),
                         legend.title=element_text(size=25),
                         legend.position=lpos,
                         strip.text.x = element_text(size = 21) )  
      pHD <- pHD + scale_color_manual(values=ourcolors) +  guides(fill=guide_legend(nrow=2,byrow=TRUE))                        
      
      
    }
    
    if(agg == "median"){
      
      # The black error bars
      # pHD <- ggplot(agg_df, aes(x=VAR, y=median, group=method)) +
      #  geom_errorbar(aes(ymin=(median-(median-first_quartile)), ymax=(median+(third_quartile-median))), width=.1) +
      #   geom_line(aes(color=method), size=3) + 
      #   #geom_point() +
      #   labs(x=xlabel, y = paste0(agg, " ", measure_title))       
      
      pHD <- ggplot(agg_df, aes(x=VAR, y=median, group=method,colour=method)) +
        #geom_errorbar(aes(ymin=Measure-se, ymax=Measure+se), width=.1,position=pd) +
        geom_errorbar(aes(ymin=(median-(median-first_quartile)), ymax=(median+(third_quartile-median))), width=.1,position=pd) +        
        geom_line(aes(color=method), size=3,position=pd) + 
        geom_point(position=pd,size=4) +
        labs(x=xlabel, y = paste0(measure_title, " (", agg, ")")) 
      pHD <- pHD + theme(axis.text.y  = element_text(size=25, colour= "black"),
                         axis.text.x  = element_text(size=25, colour= "black"),
                         axis.title.y =element_text(size=27), 
                         axis.title.x=element_text(size=27),
                         legend.text=element_text(size=22),
                         legend.title=element_text(size=25),                         
                         legend.position="top",
                         strip.text.x = element_text(size =21) )    
      pHD <- pHD + scale_color_manual(values=ourcolors)
      
    }
    
    ggsave(pHD,file=paste0(outdir,"/lineplot_", agg, "_",measure_name,".pdf"),width=13,height=10)              
  }
  
}

###############################

plot_data_barplots <- function(big_df, crash, model, measure_name,add_points) {  
  
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
    legpos="top"
    myheight=7
    xtext <- element_text(angle=90, vjust=0.5, size=14, colour= "black")
  } else if (measure_name == "Vmeasure") {
    measure_title <- "V-measure"
    fname <- "results"
    legpos="top"
    myheight=3.8  
    xtext <- element_blank()    
  } else if (measure_name == "nclusters") {
    measure_title <- "# of pred. clusters"
    fname <- "results"
    legpos="none"        
    if (one_filter) {
        myheight=3   
        xtext <- element_blank()
    } else {        
        myheight=5.5
        xtext <- element_text(angle=90, vjust=0.5, size=14, colour= "black")    
    }        
  } else if (measure_name == "nclusters2") {
    measure_title <- "# of predicted clusters"
    fname <- "results"    
    legpos="none"    
    myheight=6.3    
    xtext <- element_text(angle=90, vjust=0.5, size=14, colour= "black")    
  } else if (measure_name == "clone_prev_MAE") {
    measure_title <- "Cluster freq. MAE"
    fname <- "results"
    legpos="none"    
    #myheight=4.3
    myheight=3
    xtext <- element_blank()
  } else if (measure_name == "clone_prev_MSE") {
    measure_title <- "Clone prevalence MSE"
    fname <- "results"
    legpos="none"
    myheight=6.3 
    xtext <- element_text(angle=90, vjust=0.5, size=14, colour= "black")       
  } else if (measure_name == "uncertainty") {
    measure_title <- "Uncertainty true positive rate"
    fname <- "results"
    legpos="top"
    myheight=7    
    xtext <- element_text(angle=90, vjust=0.5, size=14, colour= "black")    
  }      
  
  xlabel = "Data set"  
      
  
    # Not used
    big_df$replicate[big_df$replicate=="0_0.95_10000_K10"] <- "filter_1"  
    big_df$replicate[big_df$replicate=="0_0.95_15000_K10"] <- "filter_2" 
    big_df$replicate[big_df$replicate=="0_0.95_20000_K10"] <- "filter_3" 
    big_df$replicate[big_df$replicate=="0_0.95_10000_K30"] <- "filter_1"  
    big_df$replicate[big_df$replicate=="0_0.95_15000_K30"] <- "filter_2" 
    big_df$replicate[big_df$replicate=="0_0.95_20000_K30"] <- "filter_3"     
    # MA: it may be 0.98
    big_df$replicate[big_df$replicate=="0_0.98_10000_K10"] <- "filter_1"  
    big_df$replicate[big_df$replicate=="0_0.98_15000_K10"] <- "filter_2" 
    big_df$replicate[big_df$replicate=="0_0.98_20000_K10"] <- "filter_3" 
    
    big_df$replicate[big_df$replicate=="0_1_0.01_K10"] <- "large_input" 
    big_df$replicate[big_df$replicate=="0_1_0.01_K30"] <- "large_input" 
  
    tmp_ind <- ((big_df$VAR == "Luo2017"))
    # tmp_ind <- ((big_df$VAR == "Luo2017")  & (big_df$replicate == "large_input"))
    print(tmp_ind)
    big_df <- big_df[!tmp_ind,]
    
    
    if (one_filter) {
        replicate_method <- paste0(as.character(big_df$method),big_df$imputed)    
    } else {        
        replicate_method <- paste0(as.character(big_df$method),"_",as.character(big_df$replicate),big_df$imputed)
    }        
    
    big_df <- cbind(big_df,replicate_method)
  
  
    big_df <- subset(big_df, replicate_method != "EpiclomalRegion_large_input *")  
    if (one_filter) {
        big_df <- subset(big_df, replicate_method != "EpiclomalRegion_filter_2")
        big_df <- subset(big_df, replicate_method != "EpiclomalRegion_filter_3")    
    }
    
    print(head(big_df))
    
    print(big_df)
    print(big_df$replicate_method)
    
#   big_df$replicate_method <- factor(big_df$replicate_method,levels= c("EpiclomalRegion_0_0.95_10000", "EpiclomalRegion_0_0.95_15000",
#                                                                         "EpiclomalRegion_0_0.95_20000" ,"EuclideanClust_0_0.95_10000" , "EuclideanClust_0_0.95_15000" ,
#                                                                         "EuclideanClust_0_0.95_20000" , "EuclideanClust_0_1_0.01"  , "DensityCut_0_0.95_10000"  ,
#                                                                         "DensityCut_0_0.95_15000"   ,   "DensityCut_0_0.95_20000"  ,    "DensityCut_0_1_0.01"  ,   "HammingClust_0_0.95_10000"  , 
#                                                                         "HammingClust_0_0.95_15000"  ,  "HammingClust_0_0.95_20000"  ,  "HammingClust_0_1_0.01"  ,      "PearsonClust_0_0.95_10000" ,   
#                                                                         "PearsonClust_0_0.95_15000"   , "PearsonClust_0_0.95_20000" ,  
#                                                                         "PearsonClust_0_1_0.01" ))
  
    # the levels are mentioned to specify the order. Don't need this though because I can relevel later
#    big_df$replicate_method <- factor(big_df$replicate_method,levels= c("EpiclomalRegion_filter_1", 
#                                                                        "EpiclomalRegion_filter_2", 
#                                                                        "EpiclomalRegion_filter_3" ,
#                                                                        "EuclideanClust_filter_1" , 
#                                                                        "EuclideanClust_filter_2", 
#                                                                        "EuclideanClust_filter_3" , 
#                                                                        "EuclideanClust_large_input", 
#                                                                        "DensityCut_filter_1"  , 
#                                                                        "DensityCut_filter_2"   ,  
#                                                                        "DensityCut_filter_3"  ,   
#                                                                        "DensityCut_large_input"  , 
#                                                                        "HammingClust_filter_1"  , 
#                                                                        "HammingClust_filter_2"  , 
#                                                                        "HammingClust_filter_3"  , 
#                                                                        "HammingClust_large_input"  , 
#                                                                        "PearsonClust_filter_1" ,   
#                                                                        "PearsonClust_filter_2"   , 
#                                                                        "PearsonClust_filter_3" ,  
#                                                                        "PearsonClust_large_input" ),
#                                                                        order = TRUE)

    big_df$replicate_method <- factor(big_df$replicate_method)   


# Do this if I use only one filter 
    if (one_filter) {
        big_df$replicate_method <- relevel(big_df$replicate_method, "EpiclomalRegion") 
    } else {
        big_df$replicate_method <- relevel(big_df$replicate_method, "EpiclomalRegion_filter_3")
        big_df$replicate_method <- relevel(big_df$replicate_method, "EpiclomalRegion_filter_2")
        big_df$replicate_method <- relevel(big_df$replicate_method, "EpiclomalRegion_filter_1")             
    }
 
    print(head(big_df))
    print(levels(big_df$VAR))
    
    Failure <- rep("F",dim(big_df)[1])
    
    big_df <- cbind(big_df,Failure)
    
    #print(big_df)   

    big_df$VAR <- factor(big_df$VAR,levels=c("Smallwood2014",
                                               "Hou2016",
#                                               "Luo2017",
                                               "Farlik2016"
#                                                "InHouse"
                                              ))
	if (measure_name != "nclusters2") {
		xlsfile <- paste0(outdir,"/SourceData.xlsx")
		write.xlsx(subset(big_df, select=c("Measure", "VAR", "replicate_method")), file=xlsfile, sheetName=measure_name, append=TRUE)	
	}
      
   if(length(grep("nclusters",measure_name)) == 0){    

    #big_df$Measure[big_df$crash == 0] <- -0.125
    big_df$Measure[big_df$crash == 0] <- 0
    big_df$Failure[big_df$crash == 1] <- ""
    #Failure[9:10] <- "F"
    
    print(big_df)
    
    #stop()
    
    
   # pHD <- ggplot(big_df, aes(x=replicate, y=Measure, fill=method)) +
   #    geom_bar(stat="identity",position=position_dodge()) +    
   #    facet_grid(~VAR) +
   #    labs(x="", y = measure_title)
    
     #geom_text(aes(label=score), size=4)
     
     
    pHD <- ggplot(big_df, aes(x=replicate_method, y=Measure, fill=method)) +
      geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.3)) + 
      geom_text(aes(label=Failure), size=4) +
      #geom_bar(width=0.4, position = position_dodge(width=0.5))
      facet_grid(~VAR,scales="free_x", space = "free_x") +
      labs(x="", y = measure_title) +
      #scale_y_continuous(breaks=c(-0.125,0.00,0.25,0.50,0.75,1.00), labels = c("Failure",0.00,0.25,0.50,0.75,1.00),limits=c(-.125,1.00)) +
      theme(plot.title = element_text(size=8), 
            axis.text.x  = xtext,
            legend.position=legpos,    
            #legend.position="top",
            axis.text.y  = element_text(size=18, colour= "black"),
            axis.title.y =element_text(size=22), 
            axis.title.x=element_text(size=20),
            legend.text=element_text(size=22) ,
            legend.title= element_text(size=22) ,
            #legend.text=element_blank() ,
            #legend.title= element_blank() ,
            strip.text.x = element_text(size = 22)
            )
                      

  pHD <- pHD + scale_fill_manual(values=ourcolors)   
             
 
  #ggsave(pHD,file=paste0(outdir,"/barplots_",measure_name,".pdf"),width=16,height=myheight)   
  ggsave(pHD,file=paste0(outdir,"/barplots_",measure_name,".pdf"),width=15,height=myheight)   
  
   }

  if(length(grep("nclusters",measure_name)) > 0) {
    big_df$Measure[big_df$crash == 0] <- -1
    
    # pHD <- ggplot(big_df, aes(x=replicate, y=Measure, fill=method)) +
    #    geom_bar(stat="identity",position=position_dodge()) +    
    #    facet_grid(~VAR) +
    #    labs(x="", y = measure_title)
    
    # MA 4Feb2019. Trying to add horizontal lines to show the correct number of clusters.
    #hline.data <- data.frame(z = c(2,2,21,6,3),  method = factor(c("Smallwood2014","Hou2016","Luo2017","Farlik2016","InHouse")))        
   
    # Camila 6Feb2019
    #hline.dat <- data.frame(VAR=levels(big_df$VAR), vl=c(2,2,21,6,3))   # For Luo
    #hline.dat <- data.frame(VAR=levels(big_df$VAR), vl=c(2,2,6,3))   # this includes Inhouse
    hline.dat <- data.frame(VAR=levels(big_df$VAR), vl=c(2,2,6))
    
    #cutoff <- data.frame( x = c(-Inf, Inf), y = 2, cutoff = factor(2) )          
    pHD <- ggplot(big_df, aes(x=replicate_method, y=Measure, fill=method)) +
      geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.3)) + 
      facet_grid(~VAR,scales="free_x", space = "free_x") +
      labs(x="", y = measure_title) +
#      scale_y_continuous(breaks=c(-1,0,2,4,6,8,10,12), labels = c("Failure",0,2,4,6,8,10,12)) +
      theme(plot.title = element_text(size=8), 
            axis.text.x  = xtext,
            #axis.text.x  = element_text(angle=90, vjust=0.5, size=14, colour= "black"),
            legend.position=legpos,
            axis.text.y  = element_text(size=18, colour= "black"),
            axis.title.y =element_text(size=22), 
            axis.title.x=element_text(size=20),
            legend.text=element_text(size=22) ,
            legend.title=element_text(size=22) ,
            strip.text.x = element_text(size = 22)
      )
    if (measure_name == "nclusters") {  
        # MA 4Feb2019: trying to add horizontal line to show the number of correct clusters
        
      pHD <- ggplot(big_df, aes(x=replicate_method, y=Measure, fill=method)) +
        geom_bar(stat="identity",width=0.5,position=position_dodge(width=0.3)) + 
     
        # Camila 6Feb2019
        geom_hline(aes(yintercept=vl), data=hline.dat,linetype="dashed", color = "black", size=1) + 
       
        facet_grid(~VAR,scales="free_x", space = "free_x") +
        labs(x="", y = measure_title) +
        #      scale_y_continuous(breaks=c(-1,0,2,4,6,8,10,12), labels = c("Failure",0,2,4,6,8,10,12)) +
        theme(plot.title = element_text(size=8), 
              axis.text.x  = xtext,
              legend.position=legpos,
              axis.text.y  = element_text(size=18, colour= "black"),
              axis.title.y =element_text(size=22), 
              axis.title.x=element_text(size=20),
              legend.text=element_text(size=22) ,
              legend.title=element_text(size=22) ,
              strip.text.x = element_text(size = 22)
        )
      
  #pHD <- pHD + geom_hline(aes(yintercept = z), hline.data, linetype="dashed", color = "coral", size=1)
   
       }        
    
    pHD <- pHD + scale_fill_manual(values=ourcolors)   
    
    
    #ggsave(pHD,file=paste0(outdir,"/barplots_",measure_name,".pdf"),width=16,height=myheight)   
    ggsave(pHD,file=paste0(outdir,"/barplots_",measure_name,".pdf"),width=15,height=myheight)   
    
  }
      
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
        #geom_errorbar(aes(ymin=Measure-se, ymax=Measure+se), width=.1,position=pd) +
        geom_errorbar(aes(ymin=(median-(median-first_quartile)), ymax=(median+(third_quartile-median))), width=.1,position=pd) +        
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



