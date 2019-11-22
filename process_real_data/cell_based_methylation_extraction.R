#======================
# libraries
#======================
# .libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos5/3.2", "/clusterapp/software/linux-x86_64-centos5/R-3.2.3/lib64/R/library"))

suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(plyr))
suppressMessages(library(IRanges))
##suppressMessages(library(reshape))
suppressMessages(library(GenomicRanges))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--data_ID", type="character", help="Some identification for the data, e.g., Smallwood2014_CGI")
parser$add_argument("--path_CpG_coordinates", type="character",help="Path to the directory containing the CpG coordinates per chromosome. Each chromosome file has to follow the format, XXX_chr1.csv")
parser$add_argument("--path_cell_data", type="character", help="Directory of the cell coverage and methylation data")
parser$add_argument("--cell_ID", type="character", help="Cell ID")
parser$add_argument("--data_type", type="character",default="bismark", help="novoalign, bismark or Farlik, default to bismark")
parser$add_argument("--genome", type="character",default="human", help="human or mouse")
parser$add_argument("--include_chrY", type="double",default="1", help="1 = include chr Y, 0 = do not include chr Y")

args <- parser$parse_args()

print(args)

#========================
# auxiliary functions
#========================

get_cov_data <- function(coverage_file, data_type) {

  cov_data <- read.table(coverage_file)

      if (data_type == "novoalign"){

        ## the columns of the bed.gz files are "chr","CpG_start","CpG_end","meth_frac","count_meth","total_reads","count_C_dimer"
        ## 6th column is total reads, so we need to transform it in count_unmeth

        cov_data <- cov_data[,1:6]
        cov_data[,6] <- cov_data[,6] -  cov_data[,5]

        cov_data[,1] <- as.character(sub("chr","",cov_data[,1]))

        #print(head(cov_data))

      }

    if (data_type == "Luo"){

      cov_data <- read.table(coverage_file,header=TRUE) ### original table already contains header

      #print(head(cov_data))
      #print(dim(cov_data))
      #print(head(CpGs_per_region))

      ##print(head(cov_data[ (cov_data[,5]/cov_data[,6]) != cov_data[,7],]))
      ##print(sum( (cov_data[,5]/cov_data[,6]) != cov_data[,7]))

      cov_data <- cov_data[!((cov_data[,5]/cov_data[,6]) != cov_data[,7]),]

      #print(dim(cov_data))

      tmp <- cov_data[,1:2]
      tmp <- cbind(tmp,cov_data[,2])
      #print(head(tmp))
      tmp <- cbind(tmp,(cov_data[,5]/cov_data[,6]),cov_data[,5],(cov_data[,6]-cov_data[,5]))
      #print(head(tmp))
      #
      cov_data <- tmp
      #
      rm(tmp)

    }

    ## Farlik data format:
    ## DNA methylation calls in the format:
    ## chromosome coordinate number-of-methylated-reads number-of-total-reads

  if (data_type == "Farlik"){

    cov_data[,1] <- sub(x=cov_data[,1],pattern="chr",replacement="")

    tmp <- cov_data[,1:2]
    tmp <- cbind(tmp,cov_data[,2])
    #print(head(tmp))
    tmp <- cbind(tmp,(cov_data[,3]/cov_data[,4]),cov_data[,3],(cov_data[,4]-cov_data[,3]))
    #print(head(tmp))

    cov_data <- tmp

    rm(tmp)

  }

  if (data_type == "scTrio") {
    cov_data[,1] <- sub(x=cov_data[,1],pattern="chr",replacement="")

    tmp <- cov_data[,1:2]
    tmp <- cbind(tmp,cov_data[,2])
    #print(head(tmp))
    tmp <- cbind(tmp,(cov_data[,8]),cov_data[,6],cov_data[,7])
    #print(head(tmp))

    cov_data <- tmp

    rm(tmp)
  }

  colnames(cov_data) <- c("chr","CpG_start","CpG_end","meth_frac","count_meth","count_unmeth")

  #print(head(cov_data))

  cov_data$chr <- as.character(cov_data$chr)

  return(cov_data)

}

extract_cov_data <- function(coverage_data,data_type,CpGs_per_region,cell_id,plot_hist,outdir){

  ### Bismark returns the coverage for both reverse and positive strands, sometimes there is coverage just for the reverse, sometimes just for the positive
  ### that's way it is important to add +1 in the end of CpG_end in our CpGs_per_region

  if(data_type == "bismark"){
    CpGs_per_region$CpG_end <- CpGs_per_region$CpG_end +  1  ### For Bismark we need to add this 1
  }

  gr1 <- makeGRangesFromDataFrame(CpGs_per_region[,1:3])
  gr2 <- makeGRangesFromDataFrame(cov_data[,1:3])
  overlap <- findOverlaps(gr1,gr2)

  #print(head(overlap))

  ### adding the coverage information to the dataframe with all CpGs in the regions of interest
  ### CpGs with no coverage will have an NA for "meth_frac","count_meth" and "count_unmeth"
  cov_columns <- matrix(NA,nrow=dim(CpGs_per_region)[1],ncol=3)
  CpGs_per_region_cov_data <- cbind(CpGs_per_region,cov_columns )

  colnames(CpGs_per_region_cov_data) <- c("chr","CpG_start","CpG_end","region_start", "region_end"  ,"region_cpgNum","region_length","region_id","meth_frac","count_meth", "count_unmeth")

  ## by doing what is commented bellow we figure out that the cov files contain data from both reverse and positive strands sometimes
  ## leading to a position in CpGs_per_region that overlap with 2 positions in the cov data
  #CpGs_per_region_cov_data[overlap@from,9:11] <- cov_data[overlap@to,4:6]
  #print(sum(!is.na(CpGs_per_region_cov_data[,9])))
  #print(sum(!is.na(CpGs_per_region_cov_data[overlap@from,9])))

  ## to address the problem above we did the following code below:
  tmp <- cov_data
  tmp$cpg_row <- NA
  tmp$cpg_row[subjectHits(overlap)] <- queryHits(overlap)
  merged <- ddply(subset(tmp, !is.na(cpg_row)), .(cpg_row), summarise, count_meth = sum(count_meth), count_unmeth = sum(count_unmeth))
  merged2 <- merged
  CpGs_per_region_cov_data[merged2$cpg_row, 9:11] <- cbind((merged2[,2]/(rowSums(merged2[,2:3]))),merged2[, 2:3])

  #print(head(CpGs_per_region_cov_data))

  ## if there is no reverse and positive strands problem then the above is the same as below"
  # tmp2 <- CpGs_per_region_cov_data
  # tmp2[queryHits(overlap),9:11] <- cov_data[subjectHits(overlap),4:6]
  # print(sum((CpGs_per_region_cov_data$meth_frac[!is.na(CpGs_per_region_cov_data$meth_frac)] != tmp2$meth_frac[!is.na(tmp2$meth_frac)])))
  # print(head(CpGs_per_region_cov_data$meth_frac[!is.na(CpGs_per_region_cov_data$meth_frac)]))
  # print(head(tmp2$meth_frac[!is.na(tmp2$meth_frac)]))

  if(plot_hist==TRUE){
    pdf(file=paste0(outdir,"hist_meth_frac_CpGs_CGI_",cell_id,".pdf"),height=7,width=7)
    hist(CpGs_per_region_cov_data$meth_frac[!is.na(CpGs_per_region_cov_data$meth_frac)],main=cell_id,xlab="Methylation fraction")
    dev.off()
  }

  sample_id <- rep(cell_id,dim(CpGs_per_region_cov_data)[1])

  CpGs_per_region_cov_data <- cbind(CpGs_per_region_cov_data,sample_id)

  colnames(CpGs_per_region_cov_data) <- c("chr","CpG_start","CpG_end","region_start", "region_end"  ,"region_cpgNum","region_length","region_id","meth_frac","count_meth", "count_unmeth","cell_id")

  return(CpGs_per_region_cov_data)

}

#====================
# END of functions
#====================

list_of_files = list.files(path = args$path_CpG_coordinates,pattern="*.tsv.gz")

print(list_of_files)

print(args$genome)

### still need to add an option when Y is not available!
if(args$genome == "human"){
  if(args$include_chrY == 1){
  chrs <- c(1:22,"X","Y")}

  if(args$include_chrY == 0){
    chrs <- c(1:22,"X")}

}

if(args$genome == "mouse"){
  chrs <- c(1:19,"X","Y")}

print(chrs)

outdir <- args$output_directory
outfile <- gzfile(paste0(outdir,"/CpG_meth_data_long_",args$data_ID,"_",args$cell_ID,".tsv.gz"))

cat(sapply(c("chr","CpG_start","CpG_end","region_start", "region_end"  ,"region_cpgNum","region_length","region_id","meth_frac","count_meth", "count_unmeth", "cell_id"), toString), file= outfile, sep="\t")
cat("\n", file= outfile, append=TRUE)

# Finding the cov file from the path and cell_id
cov_file <- Sys.glob(file.path(args$path_cell_data, paste0("*", args$cell_ID, "*")))
print("Cov file for this cell is")
print(cov_file)

if (file.exists(cov_file) && !dir.exists(cov_file))
{
  print("cov_file is a file")
  cov_data <- get_cov_data(coverage_file=cov_file, data_type=args$data_type)
}

if (file.exists(cov_file) && dir.exists(cov_file))
{
  print("cov_file is a directory")
  new_list <- list.files(cov_file)
  CpGs_per_region_cov_data_long <- NULL
}

for(c in 1:length(chrs)){

  # print(chrs[c])

  print(paste0("Chromosome ", chrs[c]))

    ## Sometimes the CpG files are called *_1.tsv.gz, sometimes they are called *_chr1.tsv.gz
    file_tmp <- list_of_files[grepl(pattern=paste0("_",chrs[c],".tsv.gz"),x=list_of_files)]
     if (length(file_tmp) == 0) {
         file_tmp <- list_of_files[grepl(pattern=paste0("_chr",chrs[c],".tsv.gz"),x=list_of_files)]
    }

    print(file_tmp)

    if (length(file_tmp) == 0 ) {
      print(paste0("There is no regions in chromosome ", chrs[c]))
    } else{

    coord_file <- paste0(args$path_CpG_coordinates,"/",file_tmp)
    # print(paste0("Coordinates file is ", coord_file))

    CpGs_per_region <- read.table(coord_file, header=TRUE) ### uncoment this please

    #print(dim(CpGs_per_region))
    #print(head(CpGs_per_region))

    if (file.exists(cov_file) && !dir.exists(cov_file))
    {
      print("cov_file is a file")

      CpGs_per_region_cov_data_per_sample <- extract_cov_data(coverage_data=cov_data,
                                                              data_type=args$data_type,CpGs_per_region=CpGs_per_region,cell_id=args$cell_ID,
                                                              plot_hist=FALSE,outdir=args$output_directory)

      #print(head(CpGs_per_region_cov_data_per_sample))
      #print(dim(CpGs_per_region_cov_data_per_sample))

      CpGs_per_region_cov_data_long <- CpGs_per_region_cov_data_per_sample

    }

    if (file.exists(cov_file) && dir.exists(cov_file))
    {

      print("cov_file is a directory")

      print(paste0(cov_file,"/",new_list[grepl(pattern=paste0("_",chrs[c],".tsv.gz"),x=new_list)]))

      cov_data <- get_cov_data(coverage_file=paste0(cov_file,"/",new_list[grepl(pattern=paste0("_",chrs[c],".tsv.gz"),x=new_list)]), data_type=args$data_type)

      CpGs_per_region_cov_data_per_sample <- extract_cov_data(coverage_data=cov_data,
                                                                  data_type=args$data_type,CpGs_per_region=CpGs_per_region,cell_id=args$cell_ID,
                                                                  plot_hist=FALSE,outdir=args$output_directory)


      CpGs_per_region_cov_data_long <- CpGs_per_region_cov_data_per_sample

      }

    #print(head(CpGs_per_region_cov_data_long))
    #print(dim(CpGs_per_region_cov_data_long))

    write.table(CpGs_per_region_cov_data_long,file= outfile,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE,quote=FALSE)

  }

  }

print("done")






