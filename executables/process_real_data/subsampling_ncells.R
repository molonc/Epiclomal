#======================
# libraries
#======================
#.libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos5/3.2", "/clusterapp/software/linux-x86_64-centos5/R-3.2.3/lib64/R/library"))

suppressMessages(library(argparse))
suppressMessages(library(pheatmap))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the directory that will contain all new input data sets")
parser$add_argument("--data_ID", type="character", help="Some identification for the data, e.g., Smallwood2014_CGI")
parser$add_argument("--file_input_Epiclomal", type="character",help="path and file to the input data for Epiclomal")
parser$add_argument("--file_final_regions", type="character",help="path and file containing the final set of regions to be used by region Epiclomal")
parser$add_argument("--file_true_membership", type="character",help="file with the true membership")
parser$add_argument("--keep_cells_prop", type="character",default="0.7_0.5_0.3",help="new number of cells for subsampling")

parser$add_argument("--number_subsamples", type="double",default="10",help="number of subsamples for each new missing proportion")
parser$add_argument("--seed", type="double",default="2",help="seed for subsampling")

args <- parser$parse_args()

print(args)

if (!is.null(args$output_directory)) {
  dir.create(args$output_directory)
}

outdir <- args$output_directory

print(outdir)

final_regions <- read.table(file=args$file_final_regions,header=TRUE)

true_membership <- read.table(file=args$file_true_membership,header=TRUE)

#print(true_membership)

input_data <- read.csv(paste0(args$file_input_Epiclomal),sep="\t",header=FALSE,na.strings=c("","NA"))

colnames(input_data) <- NULL

#print(head(input_data[1:5,1:5]))

cell_id <- as.character(input_data[,1][-1])

keep_cells_prop <- as.numeric(unlist(strsplit(args$keep_cells_prop , split ="_")))

#print(cell_id)
#print(keep_cells_prop)

set.seed(args$seed)

for(m in 1:length(keep_cells_prop)){
  
  print(m)
  print(keep_cells_prop[m])
  
  for(j in 1:args$number_subsamples){
    
    print(j)
    
    output_dir <- paste0(outdir,"/",keep_cells_prop[m],"_",j )
    
    print(output_dir)
    
    if (!dir.exists(output_dir)){
      dir.create(output_dir, showWarnings = TRUE)
    } else {
      print("This sampling already exists, going to the next one")
      next
    }     
    
    ncells_keeping <- round(length(cell_id)*keep_cells_prop[m])
    
    print(ncells_keeping)
    
    sample_ind <- sort(sample(1:length(cell_id),size=ncells_keeping,replace = FALSE))
    
    #print(sample_ind)
    #print(sort(sample_ind))
    #print(cell_id)
    print(cell_id[sort(sample_ind)])
    
    #print(dim(input_data))
    tmp <- input_data[c(1,(sample_ind + 1)),]
    #print(tmp[,1:2])
    
    #print(dim(true_membership))
    tmp_membership <- true_membership[sample_ind,]
    
    #print(tmp_membership)
    
    
    write.table(tmp, file = paste0(output_dir,"/input_Epiclomal_",args$data_ID,".tsv"), row.names = FALSE, quote = FALSE, sep = "\t", na = "")
    system(paste0("gzip --force ", output_dir,"/input_Epiclomal_",args$data_ID,".tsv"))  
    
    write.table(final_regions, file = paste0(output_dir,"/regionIDs_input_Epiclomal_",args$data_ID,".tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
    system(paste0("gzip --force ", output_dir,"/regionIDs_input_Epiclomal_",args$data_ID,".tsv"))  
    
    write.table(tmp_membership, file = paste0(output_dir,"/true_clone_membership.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
    system(paste0("gzip --force ", output_dir,"/true_clone_membership.txt"))  
    
    rm(tmp)
    rm(tmp_membership)
    
  }
  
}



print("Done!")

