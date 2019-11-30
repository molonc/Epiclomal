
.libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos5/3.2", "/clusterapp/software/linux-x86_64-centos5/R-3.2.3/lib64/R/library"))

suppressMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option

parser$add_argument("--output_dir", type="character", default="output", help="Entire output directory file")

parser$add_argument("--path_regions", type="character", help="Path to the folder containing the region files (for Farlik) or already the path to the one file with the regions")

# each replicate file should be in inputs, for example inputs/Smallwood2014_replicates.txt

args <- parser$parse_args()

print(args)

path_regions <- args$path_regions

outdir2 <- args$output_dir

if (!dir.exists(path_regions)){

  all_regions_files <- path_regions
  print(all_regions_files)

}else{

  print("path_regions is a directory")

  all_regions_files <- list.files(path = regionsFile)

  print(all_regions_files)

}

print(length(all_regions_files))

#regionsFile <- "~/Documents/shahlab15/BS-seq/whole_genome_single_cell/EPI-89/Farlik2016_regions_Fig_3f/original_regions/"
#outdir2 <- "~/Documents/shahlab15/BS-seq/whole_genome_single_cell/EPI-89/Farlik2016_regions_Fig_3f/regions_Fig_3f_Farlik2016_with_overlap/"

number_of_regions <- NULL

cat(sapply(c("chr","start","end"), toString), file= paste0(outdir2,"/all_regions_with_overlap.tsv"), sep="\t")
cat("\n", file=paste0(outdir2,"/all_regions_with_overlap.tsv") , append=TRUE)

for (j in 1:length(all_regions_files)){

  print(j)

  if (!dir.exists(path_regions)){
  region_coordinates <- read.table(all_regions_files[j],sep="\t",header=FALSE)
  }else{
    region_coordinates <- read.table(paste0(regionsFile,all_regions_files[j]),sep="\t",header=FALSE)
  }

  print(head(region_coordinates))

  if(dim(region_coordinates)[2] > 3){
    region_coordinates <- region_coordinates[,1:3]
  }

  print(head(region_coordinates))

  colnames(region_coordinates) <- c("chr","start","end")

  region_coord <- region_coordinates

  rm(region_coordinates)

  region_coord$chr <- as.character(region_coord$chr)

  rownames(region_coord) <- NULL

  print(head(region_coord))

  #print(unique(region_coord$chr))

  #region_coord <- subset(region_coord, chr %in% paste0("chr", c(1:22, "X", "Y")))

  print(unique(region_coord$chr))

  print(length(unique(region_coord$chr)))

  region_coord$chr <- factor(region_coord$chr, level = paste0("chr", c(1:22)))

  region_coord <-region_coord[order(region_coord$chr, region_coord$start), ]

  region_coord$chr <- as.character(region_coord$chr)

  list_chr <-  unique(region_coord$chr)

  print(list_chr)

  all_regions_overlap <- NULL

  for(c in 1:length(list_chr)){

    print(c)

    tmp <- region_coord[region_coord$chr == list_chr[c],]

    #print(dim(tmp))

    if(dim(tmp)[1] == 1){
      new_regions <- tmp[1,]
    }else{

    new_regions <- NULL

    current_end <- tmp$end[1]

    new_regions <- tmp[1,]

    for(i in 2:nrow(tmp)){

      print(i)
      if (tmp$start[i] < new_regions$end[dim(new_regions)[1]]){
        new_regions$end[dim(new_regions)[1]] <- tmp$end[i]
      }else{
        new_regions <- rbind(new_regions,tmp[i,])
      }


    }

    }

    all_regions_overlap <- rbind(all_regions_overlap,new_regions)

  }

  print(head(all_regions_overlap))

  ### checking overlap

  tmp <- all_regions_overlap

  tmp$overlap <- FALSE
  current_end <- tmp$end[1]
  chr_chr <- "chr1"

  for (i in 2:nrow(tmp)) {
    #print(i)
    current_start = tmp$start[i]
    if (current_end >= current_start) {
      tmp$overlap[i] <- TRUE
    } else {
      tmp$overlap[i] <- FALSE
      current_end <- tmp$end[i]
    }

    if (chr_chr == tmp$chr[i]) {

    } else {
      current_end <- tmp$end[i]
      tmp$overlap[i] <- FALSE
      chr_chr <- tmp$chr[i]
    }
  }

  print(sum(tmp$overlap==TRUE))

  number_of_regions <- c(number_of_regions,dim(all_regions_overlap)[1])

  #pdf(paste0(outdir2,all_regions_files[j],"_Histogram_length_regions_with_overlap.pdf"))
  #hist((all_regions_overlap$end-all_regions_overlap$start),main="",xlab="coord_end - coord_start")
  #dev.off()

  write.table(all_regions_overlap[,1:3], paste0(outdir2,"/all_regions_with_overlap.tsv") , sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)


}

#### Now from the file containing all regions do the following:

all_regions <- read.table(paste0(outdir2,"/all_regions_with_overlap.tsv"),sep="\t",header=TRUE)

print(head(all_regions))

print(dim(all_regions))

region_coord <- all_regions

print(length(unique(paste0(region_coord$chr,":",region_coord$star,":",region_coord$end))))

region_coord <- region_coord[!duplicated(paste0(region_coord$chr,":",region_coord$star,":",region_coord$end)),]

region_coord$chr <- factor(region_coord$chr, level = as.character(unique(region_coord$chr)) )

region_coord <-region_coord[order(region_coord$chr, region_coord$start), ]

region_coord$chr <- as.character(region_coord$chr)

list_chr <-  unique(region_coord$chr)

all_regions_overlap <- NULL

for(c in 1:length(list_chr)){

  print(c)

  tmp <- region_coord[region_coord$chr == list_chr[c],]

  print(dim(tmp))

  if(dim(tmp)[1] == 1){
    new_regions <- tmp[1,]
  }else{


    new_regions <- NULL

    current_end <- tmp$end[1]

    new_regions <- tmp[1,]

    for(i in 2:nrow(tmp)){

      #print(i)

      if (tmp$start[i] <= new_regions$end[dim(new_regions)[1]]){
        new_regions$end[dim(new_regions)[1]] <- tmp$end[i]
      }else{
        new_regions <- rbind(new_regions,tmp[i,])
      }


    }

  }

  all_regions_overlap <- rbind(all_regions_overlap,new_regions)

}

### stopped here 4:16pm May 9th

### checking overlap

tmp <- all_regions_overlap
#tmp$chr <- factor(tmp$chr, level = paste0("chr", c(1:22, "X", "Y")))
#tmp <- tmp[order(tmp$chr, tmp$start), ]
tmp$overlap <- FALSE
current_end <- tmp$end[1]
chr_chr <- "chr1"

for (i in 2:nrow(tmp)) {
  print(i)
  current_start = tmp$start[i]
  if (current_end >= current_start) {
    tmp$overlap[i] <- TRUE
  } else {
    tmp$overlap[i] <- FALSE
    current_end <- tmp$end[i]
  }

  if (chr_chr == tmp$chr[i]) {

  } else {
    current_end <- tmp$end[i]
    tmp$overlap[i] <- FALSE
    chr_chr <- tmp$chr[i]
  }
}

print("Final number of overlap: has to be zero")
print(sum(tmp$overlap==TRUE))

write.table(all_regions_overlap[,1:3], paste0(outdir2,"/all_regions_final.tsv") , sep="\t",row.names=FALSE, quote=FALSE)

dim(all_regions_overlap)
## 234700       3

#pdf(paste0(outdir2,"_Histogram_region_size_final_set.pdf"))
#hist((all_regions_overlap[,3] - all_regions_overlap[,2]) ,main="",xlab="size of regions")
#dev.off()


