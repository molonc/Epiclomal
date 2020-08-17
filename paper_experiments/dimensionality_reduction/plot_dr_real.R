suppressMessages(library(NMF))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(Rtsne))

#### I have to be here
## setwd("~/shahtemp/Epiclomal/Epiclomal/paper_experiments/dimensionality_reduction/")


num_tsne_iter <- 5000
outdir <- "../figures"

#data_id <- "Smallwood2014"
data_id <- "Hou2016"
#data_id <- "Farlik2016"

if (data_id == "Smallwood2014") {
  #input_file <- "../process_real_data/Smallwood2014/process_Smallwood2014/final_stats/0_0.95_10000/final_mean_meth_region_process_Smallwood2014.tsv"
  nmf_file <- "Smallwood2014_0_0.95_10000/nmf_rank4.rds"
  true_file <- "../process_real_data/Smallwood2014/data_Smallwood2014/true_clone_membership.txt.gz"
  epi_file <- "../process_real_data/Smallwood2014/data_Smallwood2014/true_clone_membership.txt.gz"
  #epi_file <- "../real_data/Smallwood2014/runs_epiclomal/0_0.95_10000/result_region/DIC_LINE_ELBOW_gainthr0.05_0.02/best_epi_run/cluster_MAP.tsv.gz"
  label <- paste0(data_id, "_0_0.95_10000_nmf_rank4")
  perp <- 9
} else if (data_id == "Hou2016") {
  nmf_file <- "Hou2016_0_0.95_10000/nmf_rank2.rds"
  true_file <- "../process_real_data/Hou2016/data_Hou2016/true_clone_membership.txt.gz"
  epi_file <- "../process_real_data/Hou2016/data_Hou2016/true_clone_membership.txt.gz"  
  #epi_file <- "../real_data/Hou2016/runs_epiclomal/0_0.95_10000/result_region/DIC_LINE_ELBOW_gainthr0.05_0.02/best_epi_run/cluster_MAP.tsv.gz"
  label <- paste0(data_id, "_0_0.95_10000_nmf_rank2")
  perp <- 8  
} else if (data_id == "Farlik2016") {
  nmf_file <- "Farlik2016_0_0.98_10000/nmf_rank2.rds"
  true_file <- "../process_real_data/Farlik2016/data_Farlik2016/true_clone_membership_6clusters.txt.gz"
  epi_file <- "../process_real_data/Farlik2016/data_Farlik2016/true_clone_membership_6clusters.txt.gz"
  #epi_file <- "../real_data/Farlik2016/runs_epiclomal/0_0.98_10000/result_region/DIC_LINE_ELBOW_gainthr0.05_-100/best_epi_run/cluster_MAP.tsv.gz"
  label <- paste0(data_id, "_0_0.98_10000_nmf_rank2")
  perp <- 8    
}

true <- read.csv(true_file,sep="\t",header=TRUE,check.names=FALSE)
epi <- read.csv(epi_file,sep="\t",header=TRUE,check.names=FALSE)
res <- readRDS(file=nmf_file)
htsne <- Rtsne(t(coef(res)), perplexity=perp, max_iter = num_tsne_iter, check_duplicates=FALSE)
mycols <- epi[,2]
mycols[mycols==0] <- "seagreen"
mycols[mycols==1] <- "orange3"
mycols[mycols==2] <- "royalblue"
mycols[mycols==3] <- "mediumpurple"
mycols[mycols==4] <- "steelblue"
mycols[mycols==5] <- "midnightblue"
mycols[mycols==6] <- "tomato2"


png(file=paste0(outdir,"/",label,"_tSNE_perp",perp,".png"))
par(mar=c(5.1, 4.1, 4.1, 9.1), xpd=TRUE)
plot(htsne$Y, col=mycols, pch=as.numeric(factor(true$epigenotype_id)), 
     xlab="tSNE dimension 1", ylab="tSNE dimension 2", lwd = 3, cex=1.5, 
     cex.lab=1.5, cex.axis=1.1, cex.main=1.5, cex.sub=1.5)

if (data_id == "Smallwood2014") {
  leg <- legend("topright", inset=c(-0.45,0), legend=c("2i", "SER"),
                col=c("orange3","royalblue"), 
                pch=c(1,2), lty = c(NA,NA), horiz=F, 
                lwd = 3, bty = 'n', cex=1.3)    
#  leg <- legend("topright", inset=c(-0.45,0), legend=c("2i", "SER", "Epicl. 1", "Epicl. 0"),
#                col=c("black","black","orange3", "seagreen"), 
#                pch=c(1,2,NA,NA), lty = c(NA,NA, 1, 1), horiz=F, 
#                lwd = 3, bty = 'n', cex=1.3)   
  title(main=("Smallwood2014, NMF (rank 4) + tSNE (perp. 9)"),cex.main=1.5)
} else if (data_id == "Hou2016") {
  leg <- legend("topright", inset=c(-0.45,0), legend=c("Subpop. I", "Subpop. II"),
                col=c("orange3","royalblue"), 
                pch=c(1,2), lty = c(NA,NA), horiz=F, 
                lwd = 3, bty = 'n', cex=1.3)     
#  leg <- legend("topright", inset=c(-0.45,0), legend=c("Subpop. I", "Subpop. II", "Epicl. 0", "Epicl. 1"),
#                col=c("black","black","orange3", "seagreen"), 
#                pch=c(1,2,NA,NA), lty = c(NA,NA, 1, 1), horiz=F, 
#                lwd = 3, bty = 'n', cex=1.3)   
  title(main=("Hou2016, NMF (rank 2) + tSNE (perp. 8)"),cex.main=1.5)
} else if (data_id == "Farlik2016") {
  leg <- legend("topright", inset=c(-0.45,0), legend=c("HSC", "CMP", "GMP", "MLP0", "CLP/MPP", "CLP/MLP0"),
                col=c("orange3","royalblue","mediumpurple","steelblue","midnightblue", "tomato2"), 
                pch=c(1,2,3,4,5,6),  lty = c(NA,NA,NA,NA,NA,NA), horiz=F, 
                lwd = 3, bty = 'n', cex=1.3)     
#  leg <- legend("topright", inset=c(-0.45,0), legend=c("HSC", "CMP", "GMP", "MLP0", "CLP/MPP", "CLP/MLP0", "Epicl. 2", "Epicl. 3", "Epicl. 6", "Epicl. 0", "Epicl. 5", "Epicl. 1", "Epicl. 4"),
#                col=c("black","black","black","black","black","black","orange3", "seagreen","royalblue","mediumpurple","steelblue","midnightblue","tomato2"), 
#                pch=c(1,2,3,4,5,6,NA,NA,NA,NA,NA,NA,NA), lty = c(NA,NA,NA,NA,NA,NA, 1, 1,1,1,1,1,1), horiz=F, 
#                lwd = 3, bty = 'n', cex=1.3)   
  title(main=("Farlik2016, NMF (rank 2) + tSNE (perp. 8)"),cex.main=1.5)
}

dev.off()


