suppressMessages(library(NMF))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(Rtsne))
suppressMessages(library(dbscan))
suppressMessages(library(umap))
suppressMessages(library(pheatmap))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--input_mean_meth_file", type="character", help="Input large set of mean meth by region file")
parser$add_argument("--true_membership_file", type="character", help="Input file with true methylation")
parser$add_argument("--data_ID", type="character", help="Some identification for the data, e.g., Smallwood2014_CGI")
# outdir will be the same as data_id

# global parameters

min_neighb_umap <- 10
max_neighb_umap <- 50
n_comp_umap <- c(2,10,20)
min_nmf_rank <- 2
max_nmf_rank <- 12       #15
num_nmf_runs <- 100
num_tsne_iter <- 5000
min_dbscan_points <- 5

args <- parser$parse_args()
print(args)

input_file <- args$input_mean_meth_file
true_file <- args$true_membership_file
data_id <- args$data_ID

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

calc_vmeasure <- function(cl, clusterer) {
    # make the predicted data frame
    pred <- data.frame("cell_id"=colnames(input), "epigenotype_id"=cl$cluster)
    
    # write the clusters file
    clfile <- paste0(outdir,"/",data_id,"_", clusterer ,"_clusters.csv")
    write.table(pred, file=clfile, sep="\t", col.names=TRUE, quote=FALSE,row.names=FALSE)
    system(paste0("gzip --force ", clfile))
    
    # evaluate
    rfile <- paste0(outdir,"/results_", clusterer ,".txt")
    print(paste0("Calling evaluation software for ", clusterer))
    #command <- paste0("/home/mandronescu/.local/centos6/anaconda3/bin/python /shahlab/mandronescu/MySoftware/Epiclomal/epiclomal/evaluate_clustering.py --true_clusters_file ", true_file, " --true_prevalences None --predicted_clusters_file ", clfile, ".gz --clusters_are_probabilities False --results_file ", rfile)
    command <- paste0("python ../../epiclomal/evaluate_clustering.py --true_clusters_file ", true_file, " --true_prevalences None --predicted_clusters_file ", clfile, ".gz --clusters_are_probabilities False --results_file ", rfile)
    print(command)
    system(command)
    
    vmeasure <- read.csv(rfile,sep="\t",header=TRUE,check.names=FALSE)$best_vmeasure
    return(vmeasure)
}

#######################################################
# Run UMAP
# UMAP hyper-parameters, from https://umap-learn.readthedocs.io/en/latest/parameters.html
# n_neighbors - local vs global structure - from 10 to 1/4 of the data
# min_dist - controls how tightly the points are packed together (for clustering we want a low value)
# n_components - number of dimensions to reduce to, could be more than 2 especially if we do clustering, e.g. 10 or 20
# metric
# How to use UMAP for clustering https://umap-learn.readthedocs.io/en/latest/clustering.html

compute_umap <- function(all_dat, n_neighbours, n_components) {
  custom.config = umap.defaults
  custom.config$n_neighbors = n_neighbours  # could be 30, but it has to be <= number of cells, and Hou2016 has 25 cells
  custom.config$min_dist = 0.01
  custom.config$n_components = n_components
  custom.config$random_state = 42
  
  umap(all_dat, config = custom.config)
}


for (n_comp in n_comp_umap) {
    for (n_neighb in unique(floor(seq(min_neighb_umap, ncol(input),length.out=5)))) {  
        clusterer <- paste0("UMAP_nc", n_comp, "_nn", n_neighb)
        rfile <- paste0(outdir,"/results_", clusterer ,".txt")
        if (!file.exists(rfile)) {
            print(paste0("Running UMAP with n_components ", n_comp, " and n_neighb ", n_neighb))
            dat.umap <- compute_umap(t(input), n_neighb, n_comp)
            cl <- hdbscan(dat.umap$layout, minPts = min_dbscan_points)
            vmeasure <- calc_vmeasure(cl, clusterer)
            if (n_comp == 2) {
                png(file=paste0(outdir,"/",data_id,"_UMAP_nc", n_comp, "_nn", n_neighb,".png"))          
                plot(dat.umap$layout, col=cl$cluster+1, pch=as.numeric(factor(true$epigenotype_id)))
                title(main=paste0(data_id, ", UMAP, V=", format(round(vmeasure, 2), nsmall = 2)))
                dev.off()
            }    
            print(" ...done.")
        } else {
            print(paste0("Results already exist in file ", rfile))
        }
    }        
}

#######################################################
# Now do only tSNE

# need this for tsne
max_perplexity = floor((ncol(input)-1)/3)
#max_perplexity = floor(ncol(input)-1)

# A very nice tutorial on tSNE https://distill.pub/2016/misread-tsne/
# This article explains why not to use a distance-based or density-based 
#   clustering algorithm after tsne: because the density and distances are intentionally lost
# https://stats.stackexchange.com/questions/263539/clustering-on-the-output-of-t-sne
# For tSNE, each row is an observation, each column is a variable
# Rtsne already normalizes the data by default and used some pca too
# Parameters to consider for tSNE
# - max_iter - at least 5000, should reach convergence/stability
# - perplexity (maximum is floor((ncol(input)-1)/3)) - the tradeoff between local and global
#     - shouldn't be too small (2 or 5 is too small)

for (perp in unique(floor(seq(8,max_perplexity,length.out=5)))) {
    clusterer <- paste0("tSNE_perp",perp)
    rfile <- paste0(outdir,"/results_", clusterer ,".txt")
    if (!file.exists(rfile)) {    
        print(paste0("Running only tSNE with perplexity ", perp))  
        htsne <- Rtsne(t(input), perplexity=perp, max_iter = num_tsne_iter, check_duplicates=FALSE)
        cl <- hdbscan(htsne$Y, minPts = min_dbscan_points)
        vmeasure <- calc_vmeasure(cl, clusterer)
        png(file=paste0(outdir,"/",data_id,"_tSNE_perp",perp,".png"))
        plot(htsne$Y, col=cl$cluster+1, pch=as.numeric(factor(true$epigenotype_id)))
        #title(main=paste0(data_id, ", tSNE (perp", perp, "), V=", format(round(vmeasure, 2), nsmall = 2)))
        title(main=paste0(data_id, ", tSNE (perp", perp, ")"))
        dev.off()
    } else {
        print(paste0("Results already exist in file ", rfile))
    }        
}     

#######################################################
# Run NMF on a range of ranks from 2 to 15

for (nmf_rank in seq (min_nmf_rank, max_nmf_rank)) {
    nmf_file <- paste0(outdir,"/nmf_rank", nmf_rank, ".rds")
    if (!file.exists(nmf_file)) {    
        print (paste0("Running NMF with rank ", nmf_rank))
        res <- nmf(input, nmf_rank, nrun=num_nmf_runs, .opt='vP')
        saveRDS(res, nmf_file)
    } else {
        res <- readRDS(file=nmf_file)
    }   
    
    # get all the features and plot a heatmap
    s <- extractFeatures(res)
    # do this a few times with less and more features
    for (max_num_features in c(50,100,500,1000)) {
        region_set <- c()
        for (i in 1:nmf_rank){
            if (!is.na(s[[i]])) {
                region_set <- c(region_set, s[[i]][1:min(max_num_features,length(s[[i]]))])
            }        
        }
        region_set <- unique(sort(region_set))
        selected_input <- input[region_set,]
        write.table(rownames(selected_input),file=paste0(outdir,"/NMF_regions_rank",nmf_rank,"_maxfeatures", max_num_features, ".tsv"),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)     
        pheatmap(t(selected_input),cluster_rows = TRUE,cluster_cols=FALSE,fontsize = 8, 
                 fontsize_row=2,fontsize_col=6,show_colnames = TRUE,
                 filename = paste0(outdir,"/heatmap_NMF_few_regions_mean_meth_rank", nmf_rank, "_maxfeatures", max_num_features, ".png"))               
    }    
    
    # Follow with tSNE, for up to 5 perplexities
    for (perp in unique(floor(seq(8,max_perplexity,length.out=5)))) {
        clusterer <- paste0("NMF_rank",nmf_rank,"_tSNE_perp",perp)
        rfile <- paste0(outdir,"/results_", clusterer ,".txt")
        if (!file.exists(rfile)) {            
            print(paste0("    tSNE with perplexity ", perp))
            htsne <- Rtsne(t(coef(res)), perplexity=perp, max_iter = num_tsne_iter, check_duplicates=FALSE)
            cl <- hdbscan(htsne$Y, minPts = min_dbscan_points)
            vmeasure <- calc_vmeasure(cl, clusterer)
            png(file=paste0(outdir,"/",data_id,"_NMF_rank",nmf_rank,"_tSNE_perp",perp,".png"))
            plot(htsne$Y, col=cl$cluster+1, pch=as.numeric(factor(true$epigenotype_id)))
            title(main=paste0(data_id, ", NMF (rank ",nmf_rank,")+tSNE (perp ",perp, "), V=", format(round(vmeasure, 2), nsmall = 2)))
            dev.off()
        } else {
            print(paste0("Results already exist in file ", rfile))
        }            
    }        
}    

# If I get this error, it means that we have a row full of zeros or NA
# Error: NMF::nmf - Input matrix x contains at least one null or NA-filled row.
# If some entries are NA and rest per row are not all 0, nmf works (with a warning), 
#   but the H matrix is full of NA






