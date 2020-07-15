suppressMessages(library(NMF))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(Rtsne))
suppressMessages(library(dbscan))
suppressMessages(library(umap))
suppressMessages(library(pheatmap))
library(xlsx)

min_neighb_umap <- 10
max_neighb_umap <- 50
n_comp_umap <- c(2,10,20)
min_nmf_rank <- 2
max_nmf_rank <- 12       #15
num_nmf_runs <- 100
num_tsne_iter <- 5000
min_dbscan_points <- 5

dataset <- "Smallwood2014"

if (dataset == "Smallwood2014") {
    input_file <- "../process_real_data/Smallwood2014/process_Smallwood2014/final_stats/0_0.95_10000/final_mean_meth_region_process_Smallwood2014.tsv"
    true_file <- "../process_real_data/Smallwood2014/data_Smallwood2014/true_clone_membership.txt.gz"
    data_id <- dataset
}

input_file <- "../EPI-112_inhouse_data/OUTPUT_process_INHOUSE/RUN/process_INHOUSE_process_real_data/outputs/final_stats/0_0.95_15000/final_mean_meth_region_process_INHOUSE.tsv"
true_file <- "../EPI-112_inhouse_data/data_process_INHOUSE/true_clone_membership.txt.gz"
data_id <- "InHouse"

xlsxfile <- "../EPI-109_plots_for_paper/ALL_PLOTS/SourceDataFig6.xlsx"

outdir <- data_id
dir.create(outdir)

true <- read.csv(true_file,sep="\t",header=TRUE,check.names=FALSE)

imputed_file <- paste0(outdir,"/data_imputed.csv")
if (file.exists(paste0(imputed_file,".gz"))) {
    print ("Reading the imputed file")
    input <- read.csv(paste0(imputed_file,".gz"),sep="\t",header=TRUE,check.names=FALSE)
    print (" ... done.")
} else {
    print(paste0("Reading input file ", input_file))
    input <- read.csv(input_file,sep="\t",header=TRUE,check.names=FALSE)
    # This gives features as rows and cells as columns -- this is the right input for NMF
    print(" ... done.")
    
    #input <- input[1:1000,1:200]
    
    # replace with average values, for each row
    
    print("Replacing NAs with average values")
    for (i in seq(1:nrow(input))) {
        # for some reason mean(input[i,],na.rm=TRUE) doesn't work
        vec <- input[i,!is.na(input[i,])]
        mean <- sum(vec)/length(vec)
        input[i,is.na(input[i,])] <- mean
    }
    print(" ... done.")
    
    # eliminate the empty rows (features)
    input <- input[ rowSums(input)!=0, ] 
    write.table(input, file=imputed_file, sep="\t", col.names=TRUE, quote=FALSE,row.names=TRUE)
    system(paste0("gzip --force ", imputed_file))    
}

#######################################################

# do tSNE

for (perp in c(20)) {
#for (perp in c(10, 20, 100, 150)) {
    clusterer <- paste0("tSNE_perp",perp)
    mycols <- true$epigenotype_id
    mycols[mycols==1] <- "orange3"
    mycols[mycols==2] <- "seagreen"
    mycols[mycols==3] <- "royalblue"

    tsne_file <- paste0(outdir,"/",data_id,"_tSNE_perp",perp,".rda")
    if (file.exists(tsne_file) ) {
        print("Loading tSNE data")
        load(tsne_file)
    } else {
        print(paste0("Running only tSNE with perplexity ", perp))      
        tSNE <- Rtsne(t(input), perplexity=perp, max_iter = num_tsne_iter, check_duplicates=FALSE)        
        save (tSNE, file=tsne_file)
    }    
    pdf(file=paste0(outdir,"/",data_id,"_tSNE_perp",perp,".pdf"))
    print("Printing into the excel file")
    write.xlsx(data.frame(tSNE$Y,true$epigenotype_id), file=xlsxfile, sheetName=paste0("Fig6b_", args$name), append=TRUE)
    # For some reason the line above didn't work from this script, but I ran it from R directly and it worked.
    plot(tSNE$Y, col=mycols, pch=as.numeric(factor(true$epigenotype_id)), 
         xlab="tSNE dimension 1", ylab="tSNE dimension 2", lwd = 3, cex=1.5, 
         cex.lab=1.5, cex.axis=1.1, cex.main=1.5, cex.sub=1.5)
    legend("bottomright", legend=c("SA501 TNBC", "SA532 ER+PR-Her2+", "SA609 TNBC", "Epiclomal cl. 1", "Epiclomal cl. 2", "Epiclomal cl. 3"),
           col=c("black","black","black","orange3", "seagreen", "royalblue"), 
           pch=c(1,2,3,NA,NA,NA), lty = c(NA,NA,NA, 1, 1, 1), horiz=F, 
           #inset = c(0.05, 0.05), 
           lwd = 3, bty = 'n', cex=1.3)
    title(main=paste0("tSNE visualization for the InHouse data set"),cex.main=1.5)
    dev.off()
}     

