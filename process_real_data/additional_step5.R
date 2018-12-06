#======================
# libraries
#======================
.libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos5/3.2", "/clusterapp/software/linux-x86_64-centos5/R-3.2.3/lib64/R/library"))

suppressMessages(library(argparse))
suppressMessages(library(pheatmap))
suppressMessages(library(data.table))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--data_ID", type="character", help="Some identification for the data, e.g., Smallwood2014_CGI")
parser$add_argument("--path_temp_file", type="character",help="Path to the directory containing the methylation data for each CpG in the reference genome for the regions considered. Each cell has its own tsv file")
parser$add_argument("--file_final_regions", type="character",help="File containing the final set of regions to be used to generate CpG data for Epiclomal")
parser$add_argument("--path_post_processed_CpG_data", type="character",help="Path to the directory containing the methylation data for each CpG in the reference genome for the regions considered. Each cell has its own tsv file")


args <- parser$parse_args()

if (!is.null(args$output_directory)) {
    dir.create(args$output_directory)
}

print(args)

all_CpG_cell_files <- list.files(args$path_post_processed_CpG_data)

outdir <- args$output_directory

final_regions <- read.csv(file=args$file_final_regions,sep="\t",header=FALSE)

print(head(final_regions))

colnames(final_regions) <- "region_id"

final_regions$region_id <- as.character(final_regions$region_id)

print("getting rid of CpGs with no data across all cells")

#print("loading only the data with columns with some data, that is, excluding the CpGs with no data")

tmp4 <- fread(paste0(outdir,"/input_Epiclomal_",args$data_ID,"_temp_file.tsv"),showProgress=FALSE,sep="\t",header=TRUE)

print("finding which CpGs have no data")

CpG_with_data <- rep(0,(dim(tmp4)[2]-1))

#ptm <- proc.time()

index <- apply(tmp4[,-1],2,function(x){sum(!is.na(x))})

CpG_with_data[index != 0] <- 1

print("total number of CpGs")
print(length(CpG_with_data))
print("number of CpGs with data")
print(sum(CpG_with_data == 1))

head(CpG_with_data)

print("end of finding CpGs with no data")

#print(proc.time() - ptm)

print(dim(tmp4))
print(tmp4[1:2,1:15])

no_data_index <- (which(CpG_with_data == 0)+1)
print(head(no_data_index))

tmp4 <- as.data.frame(tmp4)

print(tmp4[,no_data_index[1]])

tmp4 <- tmp4[,-no_data_index]

print(dim(tmp4))
print(tmp4[1:2,1:15])

print("saving the final data without CpGs with no data")

fwrite(tmp4, file = paste0(outdir,"/input_Epiclomal_",args$data_ID,"_only_CpGs_with_Data.tsv"), row.names = FALSE, quote = FALSE, sep = "\t", na = "",showProgress=FALSE)
system(paste0("gzip --force ", outdir,"/input_Epiclomal_",args$data_ID,"_only_CpGs_with_Data.tsv"))  

print("adjusting for right set of regions")

### Getting the right set of regions now that we eliminated some of the CpGs

  tmp <- fread(paste0(args$path_post_processed_CpG_data,"/",all_CpG_cell_files[1]),showProgress=FALSE,sep="\t",header=TRUE)

  sub_tmp <- tmp[tmp$region_id %in% final_regions$region_id,]

    print(dim(sub_tmp))

    print(length(CpG_with_data))
    
    sub_tmp <- sub_tmp[(CpG_with_data == 1),]

    print(dim(sub_tmp))
    
    col_num_id <- 1:dim(sub_tmp)[1]
    num_regions <- sum(!duplicated(sub_tmp$region_id)) 
    
    region_sizes <- diff(col_num_id[!duplicated(sub_tmp$region_id)]) ### does not include the size of last region
      
    region_sizes <- c(region_sizes, (dim(sub_tmp)[1] - sum(region_sizes)))
      
       reg_coord <- NULL
       for(r in 1:sum(!duplicated(sub_tmp$region_id))){
         if(r==1){
           reg_coord <- rbind(reg_coord,c(1, region_sizes[r]))} else{
             reg_coord <- rbind(reg_coord,c(sum(region_sizes[1:(r-1)]) + 1 , sum(region_sizes[1:r]) ) ) }
       }
      
       reg_coord <- cbind((1:num_regions),reg_coord)
       colnames(reg_coord) <- c("region_id","start","end")
      
       write.table(as.matrix(reg_coord-1 ), file = paste0(outdir,"/regionIDs_input_Epiclomal_",args$data_ID,".tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
       system(paste0("gzip --force ", outdir,"/regionIDs_input_Epiclomal_",args$data_ID,".tsv"))
    
print("End of getting input ready for Luo diamonds")   







