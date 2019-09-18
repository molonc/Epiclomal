extract_mean_meth_per_cell <- function(cell_data,region_coord){
  mean_meth <- apply(region_coord,1,function(x){mean(cell_data[x[1]:x[2]],na.rm=TRUE)})
  mean_meth[is.na(mean_meth)] <- NA
  return(mean_meth)
}

load_data <- function(input_CpG_data_file, input_regions_file){
  cached_data <- gsub(".tsv.gz", ".RDa.gz", input_CpG_data_file)
  if (file.exists(cached_data) & use_cache) {
    print("loading cached data")
    load(cached_data)
  } else {
    tmp <- read.csv(input_CpG_data_file,sep="\t",header=TRUE,check.names=FALSE)
    input_CpG_data <- as.matrix(tmp[,-1])
    rownames(input_CpG_data) <- tmp$cell_id
    rm(tmp)

    # Region coordinates
    tmp <- read.csv(input_regions_file,sep="\t",header=TRUE,check.names=FALSE)
    input_regions <- as.matrix(tmp[,-1]) + 1 ## adding 1 to match R indexing - previously coordinates were for python starting on zero
    colnames(input_regions) <- c("start","end") ## input_regions gives already the columns in input_CpG_data that correspond to which regions considered in the construction of input_CpG_data
    rownames(input_regions) <- tmp$region_id
    rm(tmp)

    mean_meth_matrix <- t(apply(input_CpG_data,1,extract_mean_meth_per_cell,region_coord=input_regions))

    save(input_CpG_data, input_regions, mean_meth_matrix, file = cached_data, compress = "gzip")
  }
  return(list("input_CpG_data" = input_CpG_data, "input_regions" = input_regions, "mean_meth_matrix" = mean_meth_matrix))
}