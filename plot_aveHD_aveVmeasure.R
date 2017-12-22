

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

#summary_table_file <- "/Users/cdesouza/Documents/shahlab15/mandronescu/EPI-73_run_synthetic_pipeline/RESULTS_NCELLS/table_NCELLS_D30_R50.txt"

summary_table <- read.table(summary_table_file,header=TRUE,na.strings="NA",sep="\t")

print(summary_table)

print(class(summary_table))

print(str(summary_table))

#stop()

outdir <- args$output_dir
#outdir <- "~/Documents/shahlab15/csouza/BS-seq/whole_genome_single_cell/EPI-73"

criterion <- args$criterion

##################
### line plots ###
##################

##adding basic_true and region_true to the plot
info_file <- paste0(outdir, "/info_", args$var, ".txt")
print(info_file)
title <- system(paste0("cat ", info_file, " | grep -v DPARPROB | grep -v REGSIZE | grep -v CONFIG | grep -v GPROB | grep -v EPROB | grep -v CPREV | grep -v ", args$var, " | perl -p -e 's/\n/ /g;'"), intern=TRUE)
title <- sub("multinomial_equal","mnon_eq",title)
print (paste0("Title is ", title))
#title <- "Camila"

print(title)

tmp <- cbind(summary_table$Avg_Vmeasure_basic_true, summary_table$Avg_Vmeasure_region_true, summary_table$Avg_Vmeasure_basic, summary_table$Avg_Vmeasure_basic_munok, summary_table$Avg_Vmeasure_region, summary_table$Avg_Vmeasure_region_munok, summary_table$Avg_Vmeasure_PBAL_Bestcut, summary_table$Avg_Vmeasure_densitycut)
print(str(tmp))
print(tmp)

#if(!is.numeric(summary_table[,1])){
#  print("ok")
#}

#stop()

pdf(paste0(outdir,"/plot_aveVmeasure_basic_VS_region_",criterion,".pdf"),height=7,width=9)

#print("SUmmaryTable is")
#print(paste0("Length of SUmmaryTable is ", length(summary_table[,1])))
#print(summary_table[,1])
#summary_table[,1] <- factor(summary_table[,1])
#print("end summary table")

if(is.numeric(summary_table[,1])){
  
matplot(summary_table[,1],tmp,lty=1,type='l', 
    ylab="Average V-measure",
    xlab=paste0(colnames(summary_table)[1]," ",criterion),
    main=title,
    cex.axis=1.2,cex.lab=1.2,xaxt="n",ylim=c(0,1), lwd=c(10,8,6,4,6,4,2,4),col=c(1,6,2,4,7,8,3,5))
# also tried to set a limit on x but didn't work xlim=c(1,length(summary_table[,1]))    
# col is color: 1=black, 2=blue, 3=green, 4=red, 6=magenta, 5=??
# trying without y limit to see if it shows better
#,ylim=c(0,1))

axis(1, summary_table[,1])
# axis(1, 1:5)
# axis(1, at=1:length(summary_table[,1]), labels=c("0.33_0.33_0.34","0.2_0.5_0.3","0.45_0.45_0.1","0.49_0.49_0.02","0.8_0.1_0.1"))   #summary_table[,1])
# axis(1, at=1:5)
legend("bottomleft",c("Basic Epiclomal True", "Region Epiclomal True", "Basic Epiclomal", "Basic munok", "Region Epiclomal","Region munok","PBALclust","densitycut"),bty="n",cex=.8,col=c(1,6,2,4,7,8,3,5),lty=c(1,1,1,1,1,1,1,1),lwd=c(10,8,6,4,6,4,2,4))
grid()
# not sure how to add text
#text(pos=2,"NCELLS=100\nNLOCI=10000")

}

if(!is.numeric(summary_table[,1])){
  
  x_tmp  <- (1:dim(tmp)[1])
  matplot(x_tmp,tmp,lty=1,type='l', 
          ylab="Average V-measure",
          xlab=paste0(colnames(summary_table)[1]," ",criterion),
          main=title,
          cex.axis=1.2,cex.lab=1.2,ylim=c(0,1),xaxt="n",lwd=c(10,8,6,4,2,4),col=c(1,6,2,4,3,5))
  # also tried to set a limit on x but didn't work xlim=c(1,length(summary_table[,1]))    
  # col is color: 1=black, 2=blue, 3=green, 4=red, 6=magenta, 5=??
  # trying without y limit to see if it shows better
  #,ylim=c(0,1))
  
  #axis(1, as.character(summary_table[,1])) ### doesn't work as well - Camila
  # axis(1, 1:5)
  axis(1, at=1:length(summary_table[,1]), labels=as.character(summary_table[,1]))
  # axis(1, at=1:5)
  legend("bottomleft",c("Basic Epiclomal True", "Region Epiclomal True", "Basic Epiclomal", "Basic munok", "Region Epiclomal","Region munok", "PBALclust","densitycut"),bty="n",cex=.8,col=c(1,6,2,4,7,8,3,5),lty=c(1,1,1,1,1,1,1,1),lwd=c(10,8,6,4,6,4,2,4))
  grid()
  
  
}

dev.off()

pdf(paste0(outdir,"/plot_aveHD_basic_VS_region_",criterion,".pdf"),height=7,width=9)
tmp <- cbind(summary_table$Avg_avgHD_basic,summary_table$Avg_avgHD_region)

  if(is.numeric(summary_table[,1])){
  matplot(summary_table[,1],tmp,lty=1,type='l',lwd=c(6,4),col=c(2,4),
    ylab="Average cell-based mean hamming distance",
    xlab=paste0(colnames(summary_table)[1]," ",criterion),
    main=title,
    cex.axis=1.2,cex.lab=1.2,xaxt="n",ylim=c(0,1))
  # trying without y limit to see if it shows better
  # ,ylim=c(0,1))
  axis(1, summary_table[,1])
  grid()
  legend("topleft",c("Basic Epiclomal","Region Epiclomal"),bty="n",col=c(2,4),lty=c(1,1),lwd=c(6,4),cex=.8)
  }


  if(!is.numeric(summary_table[,1])){
  x_tmp  <- (1:dim(tmp)[1])
  matplot(x_tmp,tmp,lty=1,type='l',lwd=c(6,4),col=c(2,4),
          ylab="Average cell-based mean hamming distance",
          xlab=paste0(colnames(summary_table)[1]," ",criterion),
          main=title,
          cex.axis=1.2,cex.lab=1.2,xaxt="n",ylim=c(0,1))

  axis(1, at=1:length(summary_table[,1]), labels=as.character(summary_table[,1]))   #summary_table[,1])
  grid()
  legend("topleft",c("Basic Epiclomal","Region Epiclomal"),bty="n",col=c(2,4),lty=c(1,1),lwd=c(6,4),cex=.8)
  }

dev.off()

##################
### box plots ####
##################

### hamming distance

number_data_sets <- as.numeric(gsub("_|D", "", str_extract(summary_table_file, "_D[0-9]+_")))

variable <- as.character(summary_table[,1])

model <- c("basic","region")

initial_path_to_each_RUN <- paste0(unlist(strsplit(summary_table_file, "/FINAL"))[1],"/RUN/D_")
print (initial_path_to_each_RUN)

print ("Boxplots for hamming distance")

method <- NULL
VAR <- NULL
hamming_distance <- NULL

counts <- c(0,0)

for(m in 1:length(model)){
  for(j in 1:length(variable)){
    print(j)
    for(i in 1:number_data_sets){
      #print(i)
      t <- try(read.table(file=paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/results_",model[m],"/",criterion,"/all_hdist_bestrun_",model[m],".tsv"),sep="\t",header=TRUE))   
      if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
        print("can't find file")
        counts[j] <- counts[j]+1 }
      else {
        hD <- read.table(file=paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/results_",model[m],"/",criterion,"/all_hdist_bestrun_",model[m],".tsv"),sep="\t",header=TRUE)
        
        hamming_distance <- c(hamming_distance,hD$mean)
        
        VAR <- c(VAR,variable[j])
        
        method <- c(method,model[m]) }
      
    }
  }
}


big_hD_df <- cbind(as.data.frame(hamming_distance),VAR,method)
colnames(big_hD_df) <- c("hD","VAR","method")

str(big_hD_df)

big_hD_df$method <- factor(big_hD_df$method,levels=c("basic","region"))
big_hD_df$VAR <- factor(big_hD_df$VAR,levels=variable)

pHD <- ggplot(big_hD_df, aes(x=method, y=hD,fill=method)) +
  geom_boxplot() + facet_grid(~VAR) +
  labs(x=" ", y = "Cell-based mean hamming distance") +
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12, colour= "black"), axis.text.y  = element_text(size=15, colour= "black"),
        #panel.background = element_rect(fill="white",colour = 'black'), 
        axis.title.y =element_text(size=15), axis.title.x=element_text(size=15),
        strip.text.x = element_text(size =12) )

ggsave(pHD,file=paste0(outdir,"/boxplot_meanHD_basic_VS_region_",criterion,".pdf"),width=13.1,height=10.6)


### V-measure
print ("Boxplots for V-measure")

method <- NULL
VAR <- NULL
Vmeasure <- NULL

counts <- c(0,0)

for(m in 1:length(model)){
  for(j in 1:length(variable)){
    print(j)
    for(i in 1:number_data_sets){
      t <- try(read.table(file=paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/results_",model[m],"/",criterion,"/all_results_bestrun_",model[m],".tsv"),sep="\t",header=TRUE))   
      if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
        print("can't find file")
        counts[j] <- counts[j]+1 }
      else {
        v <- read.table(file=paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/results_",model[m],"/",criterion,"/all_results_bestrun_",model[m],".tsv"),sep="\t",header=TRUE)
        
        Vmeasure <- c(Vmeasure,v$best_vmeasure)
        
        VAR <- c(VAR,variable[j])
        
        method <- c(method,model[m]) }
      
    }
  }
}


big_Vmeasure_df <- cbind(as.data.frame(Vmeasure),VAR,method)
colnames(big_Vmeasure_df) <- c("Vmeasure","VAR","method")

str(big_Vmeasure_df)

big_Vmeasure_df$method <- factor(big_Vmeasure_df$method,levels=c("basic","region"))
big_Vmeasure_df$VAR <- factor(big_Vmeasure_df$VAR,levels=variable)

pvmeasure <- ggplot(big_Vmeasure_df, aes(x=method, y=Vmeasure,fill=method)) +
  geom_boxplot() + facet_grid(~VAR) +
  labs(x=" ", y = "V-measure") +
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=12, colour= "black"), axis.text.y  = element_text(size=15, colour= "black"),
        #panel.background = element_rect(fill="white",colour = 'black'), 
        axis.title.y =element_text(size=15), axis.title.x=element_text(size=15),
        strip.text.x = element_text(size =12) )

ggsave(pvmeasure,file=paste0(outdir,"/boxplot_Vmeasure_basic_VS_region_",criterion,".pdf"),width=13.1,height=10.6)

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









