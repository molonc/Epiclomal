

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

args <- parser$parse_args() 

print(args)

summary_table_file <- args$summary_file

#summary_table_file <- "/Users/cdesouza/Documents/shahlab15/mandronescu/EPI-73_run_synthetic_pipeline/RESULTS_NCELLS/table_NCELLS_D30_R50.txt"

summary_table <- read.table(summary_table_file,header=TRUE,na.strings="NA",sep="\t")

# print(summary_table)

outdir <- args$output_dir
#outdir <- "~/Documents/shahlab15/csouza/BS-seq/whole_genome_single_cell/EPI-73"

##################
### line plots ###
##################

info_file <- paste0(outdir, "/info_", args$var, ".txt")
print(info_file)
title <- system(paste0("cat ", info_file, " | grep -v DPARPROB | grep -v GPROB | grep -v EPROB | grep -v CPREV | grep -v ", args$var, " | perl -p -e 's/\n/ /g;'"), intern=TRUE)
title <- sub("multinomial_equal","mnon_eq",title)
print (paste0("Title is ", title))

pdf(paste0(outdir,"/plot_aveVmeasure_basic_VS_region.pdf"),height=7,width=9)
tmp <- cbind(summary_table$Avg_Vmeasure_basic,summary_table$Avg_Vmeasure_region, summary_table$Avg_Vmeasure_PBAL_Bestcut)
print(tmp)
matplot(summary_table[,1],tmp,lty=1,type='l',lwd=c(6,4,2),col=c(2,4,3),
    ylab="Average V-measure",
    xlab=colnames(summary_table)[1],
    main=title,
    cex.axis=1.2,cex.lab=1.2,xaxt="n",ylim=c(0,1))
# col is color: 1=black, 2=blue, 3=green, 4=red, 6=magenta
# trying without y limit to see if it shows better
#,ylim=c(0,1))
axis(1, summary_table[,1])
legend("bottomright",c("Basic Epiclomal","Region Epiclomal","PBALclust"),bty="n",col=c(2,4,3),lty=c(1,1,1),lwd=c(6,4,2),cex=.8)
grid()
# not sure how to add text
#text(pos=2,"NCELLS=100\nNLOCI=10000")

dev.off()

pdf(paste0(outdir,"/plot_aveHD_basic_VS_region.pdf"),height=7,width=9)
tmp <- cbind(summary_table$Avg_avgHD_basic,summary_table$Avg_avgHD_region)
matplot(summary_table[,1],tmp,lty=1,type='l',lwd=c(6,4),col=c(2,4),
    ylab="Average cell-based mean hamming distance",
    xlab=colnames(summary_table)[1],
    main=title,
    cex.axis=1.2,cex.lab=1.2,xaxt="n",ylim=c(0,1))
# trying without y limit to see if it shows better
# ,ylim=c(0,1))
axis(1, summary_table[,1])
grid()
legend("topleft",c("Basic Epiclomal","Region Epiclomal"),bty="n",col=c(2,4),lty=c(1,1),lwd=c(6,4),cex=.8)
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
      t <- try(read.table(file=paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/results_",model[m],"/all_hdist_bestrun_",model[m],".tsv"),sep="\t",header=TRUE))   
      if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
        print("can't find file")
        counts[j] <- counts[j]+1 }
      else {
        hD <- read.table(file=paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/results_",model[m],"/all_hdist_bestrun_",model[m],".tsv"),sep="\t",header=TRUE)
        
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

ggsave(pHD,file=paste0(outdir,"/boxplot_meanHD_basic_VS_region.pdf"),width=13.1,height=10.6)


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
      t <- try(read.table(file=paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/results_",model[m],"/all_results_bestrun_",model[m],".tsv"),sep="\t",header=TRUE))   
      if("try-error" %in% class(t)) { ### could have an alternativeFunction() here
        print("can't find file")
        counts[j] <- counts[j]+1 }
      else {
        v <- read.table(file=paste0(initial_path_to_each_RUN,colnames(summary_table)[1],"_",variable[j],"_",i,"_epiclomal_synthetic/outputs/results_",model[m],"/all_results_bestrun_",model[m],".tsv"),sep="\t",header=TRUE)
        
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

ggsave(pvmeasure,file=paste0(outdir,"/boxplot_Vmeasure_basic_VS_region.pdf"),width=13.1,height=10.6)



