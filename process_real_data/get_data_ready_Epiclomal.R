#======================
# libraries
#======================
# .libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos5/3.2", "/clusterapp/software/linux-x86_64-centos5/R-3.2.3/lib64/R/library"))

suppressMessages(library(argparse))
suppressMessages(library(pheatmap))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--data_ID", type="character", help="Some identification for the data, e.g., Smallwood2014_CGI")
parser$add_argument("--path_post_processed_CpG_data", type="character",help="Path to the directory containing the methylation data for each CpG in the reference genome for the regions considered. Each cell has its own tsv file")
parser$add_argument("--file_final_regions", type="character",help="File containing the final set of regions to be used to generate CpG data for Epiclomal")
parser$add_argument("--filter_CpG_no_data", type="double",default="0",help="0 = CpGs with no data will be filtered out, 1 = they are kept")

parser$add_argument("--LuoDiamond", type="double",default="0",help="0 = usual step 5, 1 = special step 5 for producing the diamonds for Luo data set")

args <- parser$parse_args()

if (!is.null(args$output_directory)) {
  dir.create(args$output_directory, showWarnings = FALSE, recursive=TRUE)
}

print(args)

outdir <- args$output_directory

all_CpG_cell_files <- list.files(args$path_post_processed_CpG_data, pattern = "*.tsv*")

# print(all_CpG_cell_files)

###########################
### Auxiliary Functions ###
###########################

binary_function <- function(x){ ### x is a vector with the methylation fraction of a given CpG across all cells
  x[x < 0.5]  <- 0
  x[x > 0.5]  <- 1
  x[x == 0.5] <- NA
  return(x)
}


### End of auxiliary functions

if(args$LuoDiamond == 1){

  #ptm <- proc.time()

  print("getting input ready for Luo diamonds")

  final_regions <- read.csv(file=args$file_final_regions,sep="\t",header=FALSE)

  print(head(final_regions))

  colnames(final_regions) <- "region_id"

  final_regions$region_id <- as.character(final_regions$region_id)

  #system.time(tmp <- read.csv(file.path(args$path_post_processed_CpG_data, all_CpG_cell_files[1]),sep="\t",header=TRUE))

  suppressMessages(library(data.table))

  system.time(tmp <- fread(file.path(args$path_post_processed_CpG_data, all_CpG_cell_files[1]),showProgress=FALSE,sep="\t",header=TRUE))

  sub_tmp <- tmp[tmp$region_id %in% final_regions$region_id,]

  id <- with(sub_tmp, paste0(chr, ":", CpG_start))

  if(args$filter_CpG_no_data == 1){

    col_num_id <- 1:dim(sub_tmp)[1]
    num_regions <- sum(!duplicated(sub_tmp$region_id))

    region_sizes <- diff(col_num_id[!duplicated(sub_tmp$region_id)]) ### does not include the size of last region
    region_sizes <- c(region_sizes, (dim(sub_tmp)[1] - sum(region_sizes)))

    reg_coord <- NULL
    for(r in 1:sum(!duplicated(sub_tmp$region_id))){
      if(r==1){
        reg_coord <- rbind(reg_coord,c(1, region_sizes[r]))} else{
        reg_coord <- rbind(reg_coord,c(sum(region_sizes[1:(r-1)]) + 1 , sum(region_sizes[1:r]) ) )
      }
    }

    reg_coord <- cbind((1:num_regions),reg_coord)
    colnames(reg_coord) <- c("region_id","start","end")

    filename = gzfile(file.path(outdir, paste0("regionIDs_input_Epiclomal_",args$data_ID,".tsv.gz")))
    write.table(as.matrix(reg_coord-1 ), file = filename, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  print(dim(sub_tmp))

  epiclomal_input_file = file.path(outdir, paste0("input_Epiclomal_",args$data_ID,".tsv"))
  cat(sapply(c("cell_id",id), toString), file = epiclomal_input_file , sep="\t")
  cat("\n", file=epiclomal_input_file, append=TRUE)

  rm(tmp)
  rm(sub_tmp)

  CpG_with_data <- rep(0,length(id))

  for (c in 1:length(all_CpG_cell_files)){
    print(c)

    ##if(c==3) break

    ### saving a file with regions for Epiclomal

    #tmp <- read.csv(file.path(args$path_post_processed_CpG_data, all_CpG_cell_files[c]),sep="\t",header=TRUE)
    tmp <- fread(file.path(args$path_post_processed_CpG_data, all_CpG_cell_files[c]),showProgress=FALSE,sep="\t",header=TRUE)
    print("  ...read")

    sub_tmp <- tmp[tmp$region_id %in% final_regions$region_id,]  # fast
    print("  ...got sub_tmp")

    cell_id <- as.character(sub_tmp$cell_id[1])   # fast
    print("  ...got cell_id")

    CpG_data <- as.vector(binary_function(x=as.matrix((sub_tmp$meth_frac))))
    print("  ...got CpG_data")

    CpG_with_data[!is.na(CpG_data)] <- 1
    print("  ...got CpG_with_data")

    #print(head(CpG_data))
    #print(head(t(as.matrix(c(cell_id,CpG_data)))))
    write.table(t(as.matrix(c(cell_id,CpG_data))), file = epiclomal_input_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t", na = "",append=TRUE)
    print("  ...wrote in table")
  }

  if(args$filter_CpG_no_data == 1){
    print("End of getting input ready for Luo diamonds")
  }

  if(args$filter_CpG_no_data == 0){

    print("getting rid of CpGs with no data across all cells")

    print("total number of CpGs")
    print(length(CpG_with_data))
    print("number of CpGs with data across all cells:")
    print(sum(CpG_with_data == 1))

    #print(head(which(CpG_with_data == 1)))
    #print(head((which(CpG_with_data == 1)+1)))

    print("loading only the data with columns with some data, that is, excluding the CpGs with no data")
    tmp4 <- fread(epiclomal_input_file,showProgress=FALSE,select=c(1,(which(CpG_with_data == 1)+1)),sep="\t",header=TRUE)

    print(dim(tmp4))
    print(tmp4[1:2,1:15])
    #print(sum(CpG_data_filter == tmp4))

    print("saving the final data without CpGs with no data")

    file.rename(from=file.path(outdir, paste0("input_Epiclomal_",args$data_ID,".tsv")), to=file.path(outdir, paste0("input_Epiclomal_",args$data_ID,"_temp_file.tsv")))
    ## file.remove(file.path(outdir, paste0("input_Epiclomal_",args$data_ID,".tsv")))

    #tmp4 <- as.data.frame(tmp4)
    print(class(tmp4))

    epiclomal_input_file = gzfile(file.path(outdir, paste0("input_Epiclomal_",args$data_ID,".tsv.gz")))
    fwrite(tmp4, file = epiclomal_input_file, row.names = FALSE, quote = FALSE, sep = "\t", na = "",showProgress=FALSE)

    print("adjusting for right set of regions")

    ### Getting the right set of regions now that we eliminated some of the CpGs

    print(dim(sub_tmp))

    print(length(CpG_with_data))

    sub_tmp <- sub_tmp[(CpG_with_data == 1),]

    print(dim(sub_tmp))

    col_num_id <- 1:dim(sub_tmp)[1]
    num_regions <- sum(!duplicated(sub_tmp$region_id))

    region_sizes <- diff(col_num_id[!duplicated(sub_tmp$region_id)]) ### does not include the size of last region

    region_sizes <- c(region_sizes, (dim(sub_tmp)[1] - sum(region_sizes)))

    reg_coord <- NULL
    for(r in 1:sum(!duplicated(sub_tmp$region_id))) {
      if(r==1){
        reg_coord <- rbind(reg_coord,c(1, region_sizes[r]))
      } else {
        reg_coord <- rbind(reg_coord,c(sum(region_sizes[1:(r-1)]) + 1 , sum(region_sizes[1:r])))
      }
    }

    reg_coord <- cbind((1:num_regions),reg_coord)
    colnames(reg_coord) <- c("region_id","start","end")

    region_file = gzfile(file.path(outdir, paste0("regionIDs_input_Epiclomal_",args$data_ID,".tsv.gz")))
    write.table(as.matrix(reg_coord-1 ), file = region_file, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  print("End of getting input ready for Luo diamonds")

  #finaltime <- proc.time() - ptm
}

if (args$LuoDiamond == 0) {
  print("getting input ready")

  #ptm = proc.time()

  final_regions <- read.csv(file=args$file_final_regions,sep="\t",header=FALSE)

  colnames(final_regions) <- "region_id"

  final_regions$region_id <- as.character(final_regions$region_id)

  print(head(final_regions))

  #print(str(final_regions))

  # some setup with first cell
  tmp <- read.csv(file.path(args$path_post_processed_CpG_data,all_CpG_cell_files[1]),sep="\t",header=TRUE)
  sub_tmp <- tmp[tmp$region_id %in% final_regions$region_id,]

  ### saving a file with regions for Epiclomal
  col_num_id <- 1:dim(sub_tmp)[1]
  num_regions <- sum(!duplicated(sub_tmp$region_id))

  region_sizes <- diff(col_num_id[!duplicated(sub_tmp$region_id)]) ### does not include the size of last region

  region_sizes <- c(region_sizes, (dim(sub_tmp)[1] - sum(region_sizes)))

  reg_coord <- NULL
  for(r in 1:sum(!duplicated(sub_tmp$region_id))){
    if(r==1){
      reg_coord <- rbind(reg_coord,c(1, region_sizes[r]))
    } else {
      reg_coord <- rbind(reg_coord,c(sum(region_sizes[1:(r-1)]) + 1 , sum(region_sizes[1:r]) ) )
    }
  }

  reg_coord <- cbind((1:num_regions),reg_coord)
  colnames(reg_coord) <- c("region_id","start","end")

  if(args$filter_CpG_no_data == 1){
    region_file = gzfile(file.path(outdir, paste0("regionIDs_input_Epiclomal_",args$data_ID,".tsv.gz")))
    write.table(as.matrix(reg_coord-1 ), file = region_file, row.names = FALSE, quote = FALSE, sep = "\t")

  }

  id <- with(sub_tmp, paste0(chr, ":", CpG_start))

  number_cells <- length(all_CpG_cell_files)
  total_number_regions_filtered <- num_regions
  number_regions_single_CpG <- sum(region_sizes == 1)

  cell_id <- character(number_cells)
  CpG_data <- matrix(, nrow = length(id), ncol = number_cells)
  mono_meth_prop <- numeric(number_cells)

  start <- 1

  cached_data <- file.path(outdir, "epiclomal_input_Luo_Diamond_0.Rda.gz")
  if (file.exists(cached_data)) {
    print("loading cached epiclomal input data")
    load(cached_data)
  }
  if (start != FALSE) {
    print(paste('starting from cell', start))

    for (c in start:number_cells) {
      tmp <- read.csv(file.path(args$path_post_processed_CpG_data, all_CpG_cell_files[c]),sep="\t",header=TRUE)

      sub_tmp <- tmp[tmp$region_id %in% final_regions$region_id,]

      #print(unique(sub_tmp$region_id) == final_regions$region_id) ### it is working!

      print(as.character(sub_tmp$cell_id[1]))

      mono_meth <- sum( (sub_tmp$meth_frac[!is.na(sub_tmp$meth_frac)] != 0) & (sub_tmp$meth_frac[!is.na(sub_tmp$meth_frac)] != 1) )

      mono_meth_prop[c] <- mono_meth/length(sub_tmp$meth_frac[!is.na(sub_tmp$meth_frac)])

      CpG_data[,c] <- binary_function(x=as.matrix((sub_tmp$meth_frac)))

      cell_id[c] <- as.character(sub_tmp$cell_id[1])

      start <- start + 1

      save(mono_meth_prop, CpG_data, cell_id, start, file = cached_data, compress = "gzip")
    }

    start <- FALSE
    save(mono_meth_prop, CpG_data, cell_id, start, file = cached_data, compress = "gzip")

  } else { print("data already processed, loaded from cache") }

  CpG_data <- t(CpG_data)
  print(dim(CpG_data))
  print(length(id))

  colnames(CpG_data) <- id

  miss_prop_per_cell <- apply(CpG_data,1,function(x){sum(is.na(x))})/(dim(CpG_data)[2])
  ave_miss_prop <- mean(miss_prop_per_cell)
  print(ave_miss_prop)

  ave_mono_meth_prop <- mean(mono_meth_prop, na.rm=TRUE)

  print(ave_mono_meth_prop)

  ### making a plot for the proportion of CpGs with methylation as a fraction amongst all CpGs with data
  pdf(file.path(outdir, paste0("hist_mono_meth_prop_filtered_",args$data_ID,".pdf")))
  hist(mono_meth_prop,main="Proportion of CpGs with methylation as a fraction",xlab="Proportion of CpGs")
  dev.off()

  CpG_with_data <- apply(CpG_data,2,function(x){sum(!is.na(x))})  ### checking how many cells have data for each CpG

  ### make a plot for this CpG_with_data
  pdf(file.path(outdir, paste0("hist_number_cells_with_data_per_CpG_filtered_",args$data_ID,".pdf")))
  hist(CpG_with_data,main="Number of cells with data across CpGs",xlab="Number of cells")
  dev.off()

  num_CpG_no_data <- sum(CpG_with_data == 0) ### checking how many CpGs have no data across all cells

  IQR_CpG_data <- apply(CpG_data,2,function(x){IQR(x,na.rm=TRUE)})

  pdf(file.path(outdir, paste0("hist_CpG_based_IQR_filtered_",args$data_ID,".pdf")))
  hist(IQR_CpG_data,main="Methylation IQR across cells per CpG",xlab="Methylation IQR")
  dev.off()

  #print(dim(CpG_data) )
  #print(length(CpG_with_data == 0))
  #print(sum(CpG_with_data == 0))


  if(args$filter_CpG_no_data == 0){

    print("getting rid of CpGs with no data across all cells")

    CpG_data <- CpG_data[,!(CpG_with_data == 0)]

    ### getting the right regions for Epiclomal after getting rid of CpGs with no data

    tmp <- read.csv(file.path(args$path_post_processed_CpG_data, all_CpG_cell_files[1]),sep="\t",header=TRUE)

    sub_tmp <- tmp[tmp$region_id %in% final_regions$region_id,]

    sub_tmp <- sub_tmp[!(CpG_with_data == 0),]

    col_num_id <- 1:dim(sub_tmp)[1]
    num_regions <- sum(!duplicated(sub_tmp$region_id))

    region_sizes <- diff(col_num_id[!duplicated(sub_tmp$region_id)]) ### does not include the size of last region

    region_sizes <- c(region_sizes, (dim(sub_tmp)[1] - sum(region_sizes)))

    reg_coord <- NULL
    for(r in 1:sum(!duplicated(sub_tmp$region_id))){
      if(r==1){
        reg_coord <- rbind(reg_coord,c(1, region_sizes[r]))
      } else {
        reg_coord <- rbind(reg_coord,c(sum(region_sizes[1:(r-1)]) + 1 , sum(region_sizes[1:r]) ) )
      }
    }

    reg_coord <- cbind((1:num_regions),reg_coord)
    colnames(reg_coord) <- c("region_id","start","end")

    filename = gzfile(file.path(outdir, paste0("regionIDs_input_Epiclomal_",args$data_ID,".tsv.gz")))
    write.table(as.matrix(reg_coord-1 ), file = filename, row.names = FALSE, quote = FALSE, sep = "\t")


  }

  print(dim(CpG_data))

  CpG_with_data <- apply(CpG_data,2,function(x){sum(!is.na(x))})
  print(sum(!(CpG_with_data == 0)))

  final_miss_prop <- apply(CpG_data,1,function(x){sum(is.na(x))/length(x)})

  print(final_miss_prop)

  print(which(final_miss_prop == 1))

  filename = gzfile(file.path(outdir, paste0("final_miss_prop_per_cell_",args$data_ID,".tsv.gz")))
  write.table(cbind(cell_id,final_miss_prop), file = filename, row.names = FALSE, col.names=FALSE, quote = FALSE, sep = "\t")

  if( (dim(CpG_data)[2]) > 250){
    pheatmap(CpG_data[,1:250],cluster_rows = FALSE,cluster_cols=FALSE, cellwidth = 5,
      cellheight = 5,fontsize = 8,
      #main = paste0("CpG-based methylation data for ", args$name),
      #gaps_row = index_gaps,fontsize_row=6,fontsize_col=4,
      #annotation_names_row = FALSE, annotation_names_col= FALSE,
      #gaps_col=(input_regions[,2][1:(R-1)] + 1),
      show_colnames=FALSE,
      #annotation_col=annotation_col,
      filename = file.path(outdir, paste0("final_sample_CpG_based_PLOT_",args$data_ID,".pdf")))
  } else {
    pheatmap(CpG_data,cluster_rows = FALSE,cluster_cols=FALSE, cellwidth = 5,
    cellheight = 5,fontsize = 8,
    #main = paste0("CpG-based methylation data for ", args$name),
    #gaps_row = index_gaps,fontsize_row=6,fontsize_col=4,
    #annotation_names_row = FALSE, annotation_names_col= FALSE,
    #gaps_col=(input_regions[,2][1:(R-1)] + 1),
    show_colnames=FALSE,
    #annotation_col=annotation_col,
    filename = file.path(outdir, paste0("final_sample_CpG_based_PLOT_",args$data_ID,".pdf")))

  }

  CpG_data <- cbind(cell_id,CpG_data)

  cat("\n", file=file.path(outdir, paste0("cells_no_data_",args$data_ID,".tsv")), append=FALSE)

  if( sum(final_miss_prop == 1) > 0 ){
    print("Warning: there are cells with no data across all CpGs")

    CpG_data <- CpG_data[-which(final_miss_prop == 1),]

    print(cell_id)
    print(which(final_miss_prop == 1))

    cells_no_data <- cell_id[which(final_miss_prop == 1)]

    print(cells_no_data)

    write.table(cells_no_data, file=file.path(outdir, paste0("cells_no_data_",args$data_ID,".tsv")) , sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

  }


  print(head(CpG_data)[,head(colnames(CpG_data))])

  epiclomal_input_file = gzfile(file.path(outdir, paste0("input_Epiclomal_",args$data_ID,".tsv.gz")))

  write.table(CpG_data, file = epiclomal_input_file, row.names = FALSE, quote = FALSE, sep = "\t", na = "")

  final_number_loci <- (dim(CpG_data)[2]-1) ### number of columns in CpG_data minus the columns corresponding to cell_id

  final_number_cells <- (dim(CpG_data)[1])

  ave_miss_prop <- mean(final_miss_prop[which(final_miss_prop != 1)])

  info_filtered <- as.matrix(c(final_number_loci,final_number_cells,total_number_regions_filtered,number_regions_single_CpG,ave_miss_prop,num_CpG_no_data,sum(final_miss_prop == 1),ave_mono_meth_prop))

  rownames(info_filtered) <- c("final number of loci","final number of cells","number of regions","number of regions containing only 1 CpG", "final average missing proportion","number of CpGs with no data","number of cells with no data","average proportion of methylation as fraction")

  print("Table with some info for filtered data")
  print(info_filtered)

  write.table(info_filtered,file.path(outdir, paste0("filtered_data_info_",args$data_ID,".tsv")),sep="\t",quote=FALSE,col.names=FALSE,append=FALSE)

  #finaltime <- proc.time() - ptm
}

print("Done")





