
#======================
# libraries
#======================
suppressMessages(library(argparse))

suppressMessages(library(data.table))

RSCRIPT <- "Rscript"

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--name_dataset", type="character", help="name of the dataset")
parser$add_argument("--input_directory", default="None", type="character", help="path to directory containing CpG data")
parser$add_argument("--output_directory", type="character", help="Path to the output directory")


args <- parser$parse_args()
print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)

all_CpG_cell_files <- list.files(args$input_directory, pattern = "*.tsv*")

print(all_CpG_cell_files)

percentage_monoallelic_all <- NULL
percentage_monoallelic_more1read <- NULL
percentage_CpGs_more1read <- NULL
percentage_CpGs_with_data <- NULL

for (c in 1:length(all_CpG_cell_files)){
  print(c)
  
  #if(c==5) break
  
  tmp <- fread(file.path(args$input_directory, all_CpG_cell_files[c]),select=c("meth_frac", "count_meth","count_unmeth"),showProgress=FALSE,sep="\t",header=TRUE)
  print("  ...read")
  #print(head(tmp))
  #print(dim(tmp))
  
  #print(unique(tmp$meth_frac))
  
  #print(sum(!is.na(tmp$meth_frac))/dim(tmp)[1])
  
  percentage_CpGs_with_data  <- c(percentage_CpGs_with_data,(sum(!is.na(tmp$meth_frac))/dim(tmp)[1])*100)
  
  
  #print(head(tmp[!is.na(tmp$meth_frac),],100))
  
  #print(dim(tmp[!is.na(tmp$meth_frac),]))
  
  print("percentage of CpGs with partial methylation among all of those with data")
  p_all <- (sum(tmp$meth_frac[!is.na(tmp$meth_frac)] > 0 & tmp$meth_frac[!is.na(tmp$meth_frac)] < 1)/dim(tmp[!is.na(tmp$meth_frac),])[1])*100
  print(p_all)
  
  #print(sum(tmp$meth_frac[!is.na(tmp$meth_frac)] > 0 & tmp$meth_frac[!is.na(tmp$meth_frac)] < 1))
  #print(sum(tmp$meth_frac[!is.na(tmp$meth_frac)] == 0))
  #print(sum(tmp$meth_frac[!is.na(tmp$meth_frac)] == 1))
  
  percentage_monoallelic_all <- c(percentage_monoallelic_all , p_all)

  print("percentage of CpGs with partial methylation among those with more than one count/read")
  #counts <- rowSums(tmp[,2:3])
  #print(head(counts,10))
  #print(head(tmp,10))
  #print(sum(!is.na(counts)))
  #print(sum(counts[!is.na(counts)] > 1))
  
  tmp2 <- tmp[!is.na(tmp$meth_frac),]
  #print(length(tmp2$meth_frac[which((tmp2$count_meth + tmp2$count_unmeth) > 1)]))
  percentage_CpGs_more1read <- c(percentage_CpGs_more1read,(length(tmp2$meth_frac[which((tmp2$count_meth + tmp2$count_unmeth) > 1)])/dim(tmp2)[1])*100)

  #print( ( sum(tmp2$meth_frac > 0 & tmp2$meth_frac < 1) ))
  p_more1read <- (sum(tmp2$meth_frac > 0 & tmp2$meth_frac < 1) / length(tmp2$meth_frac[which((tmp2$count_meth + tmp2$count_unmeth) > 1)]))*100
  print(p_more1read)
  
  percentage_monoallelic_more1read <- c(percentage_monoallelic_more1read,p_more1read)
  
  rm(tmp)
  rm(tmp2)
}

print(percentage_monoallelic_more1read )
print(percentage_monoallelic_all)

df <- data.frame(percentage_monoallelic_all, percentage_monoallelic_more1read, percentage_CpGs_more1read,percentage_CpGs_with_data)
write.csv(df, file=paste0(outdir,"/monoallelic_info_",args$name_dataset,".csv"),row.names = FALSE)

print("done")
