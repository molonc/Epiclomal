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
parser$add_argument("--keep_cells_prop", type="character",default="0.7_0.5_0.3",help="new number of cells for subsampling")

parser$add_argument("--read_size", type="character",default="10_2",help="read size")

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

input_data <- read.csv(paste0(args$file_input_Epiclomal),sep="\t",header=FALSE,na.strings=c("","NA"))
colnames(input_data) <- NULL

cell_id <- as.character(input_data[,1][-1])

input_data <- as.matrix(input_data[,-1])

CpG_coord <- as.vector(input_data[1,])

input_data <- input_data[-1,]
input_data <- apply(input_data,2,function(x){as.numeric(x)})

 print(length(cell_id))
 print(length(CpG_coord))
 print(dim(input_data))

rownames(input_data) <- cell_id
colnames(input_data) <- CpG_coord

input_data <- t(input_data)

print(dim(input_data))
head(input_data)

indices_with_data <- as.matrix(which(!is.na(input_data),arr.ind=TRUE))

#print(head(indices_with_data))
#print(tail(indices_with_data))
#print(dim(indices_with_data))

current_missing <- sum(is.na(input_data))/(dim(input_data)[1]*dim(input_data)[2])
print("Original input data missing proportion:")
print(current_missing)

new_missing <- as.numeric(unlist(strsplit(args$new_miss_prop, split ="_")))

read_size <- as.numeric(unlist(strsplit(args$read_size, split ="_")))

mean_read_size <- read_size[1]
sd_read_size <- read_size[2]

n_loci <- dim(input_data)[1]

set.seed(args$seed)

for(m in 1:length(new_missing)){
  
  print(m)
  print(new_missing[m])

  if(current_missing > new_missing[m]){
  print("Stop - original missing proportion is already larger than the one requested for subsampling")  
  }else{
  
    for(j in 1:args$number_subsamples){
      
    print(j)
  
    number_extra_missing <- round(((dim(input_data)[1]*dim(input_data)[2])*new_missing[m]) - ((dim(input_data)[1]*dim(input_data)[2])*current_missing))
  
    #print(number_extra_missing)
    
    tmp <- input_data
    
    new_extra_missing <- number_extra_missing
    
    sample_previous <- NULL

    s <- 0
    
    print("Starting while loop")
    
    while ( ( new_extra_missing >= 0 )  ){

      s <- s+1
      
      #print(s)
      #print(new_extra_missing)
      
      read_length <- round(rnorm(1,mean=mean_read_size,sd=sd_read_size))
      
      sample_ind <- sample(1:dim(indices_with_data)[1],size=1)
      
      sample_previous <- c(sample_previous,sample_ind)
      
      if(sum(sample_ind == sample_previous[-length(sample_previous)]) == 1){
        #print("Sampling another index because the current one is equal to a previous one")
        
        #i=0
        while(sum(sample_ind == sample_previous[-length(sample_previous)]) == 1){
          #i <- i+1
          #print(i)
          sample_ind <- sample(1:dim(indices_with_data)[1],size=1)
          sample_previous[length(sample_previous)] <- sample_ind
          
        }
        
      }
    
      
      new_miss_indices <- indices_with_data[sample_ind,]
      
      obs_indices <- new_miss_indices
         
         if( new_extra_missing > read_length ){
          
           mid <- round(read_length/2) 
                      
           if((obs_indices[1] - mid) >= 1){
             
             if( (n_loci - obs_indices[1]) >= ((read_length-mid)-1) ){
             tmp[((obs_indices[1]-mid):(obs_indices[1]+((read_length-mid)-1))),obs_indices[2]] <- NA
             
             }else{
             tmp[((obs_indices[1]-mid):n_loci),obs_indices[2]] <- NA   
             }
      
           }else{
             tmp[(((obs_indices[1] - mid)-(obs_indices[1] - mid)+1):((obs_indices[1])+((read_length-mid)-1))),obs_indices[2]] <- NA   
           }
           
           
         }else{
           
           #print("Read length is bigger than amount of missing data needed")
           
            if ( (n_loci - obs_indices[1]) >= ( new_extra_missing - 1 ) ) {
            tmp[(obs_indices[1]: ( obs_indices[1]+ (new_extra_missing-1))),obs_indices[2]] <- NA
              }else { 
               tmp[(obs_indices[1]:n_loci),obs_indices[2]] <- NA
              }
          
            }
               
                c_missing <- sum(is.na(tmp))/(dim(tmp)[1]*dim(tmp)[2])
                
                new_extra_missing <- round(((dim(input_data)[1]*dim(input_data)[2])*new_missing[m]) - ((dim(input_data)[1]*dim(input_data)[2])*c_missing))
                

             }


    new_data_more_missing <- t(tmp)
    rm(tmp)
    
    print("missing proportion after adding more missing data:")
    print(sum(is.na(new_data_more_missing))/(dim(new_data_more_missing)[1]*dim(new_data_more_missing)[2]))
  
    #print(length(cell_id))
    #print(dim(new_data_more_missing))
    
    new_data_more_missing <- cbind(cell_id,new_data_more_missing)
  
    colnames(new_data_more_missing)[1] <- "cell_id"  
  
    #print(head(new_data_more_missing))

    output_dir <- paste0(outdir,"/",new_missing[m],"_",j )

    #print("output directory is")
    #print(output_dir)

    if (!dir.exists(output_dir)){
    dir.create(output_dir, showWarnings = TRUE)
    }

    write.table(new_data_more_missing, file = paste0(output_dir,"/input_Epiclomal_",args$data_ID,".tsv"), row.names = FALSE, quote = FALSE, sep = "\t", na = "")
    system(paste0("gzip --force ", output_dir,"/input_Epiclomal_",args$data_ID,".tsv"))  

    write.table(final_regions, file = paste0(output_dir,"/regionIDs_input_Epiclomal_",args$data_ID,".tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
    system(paste0("gzip --force ", output_dir,"/regionIDs_input_Epiclomal_",args$data_ID,".tsv"))  

    }

}

}

print("Done!")

