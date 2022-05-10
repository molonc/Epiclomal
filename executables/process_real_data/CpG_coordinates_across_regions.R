#' ---
#' author: "Camila de Souza"
#' output:
#'   html_document:
#'     keep_md: TRUE
#' ---

### This code returns the CpG coordinates within each region in a list

#======================
# libraries
#======================
#.libPaths(c("/extscratch/shahlab/dalai/R/x86_64-pc-linux-gnu-library-centos5/3.2", "/clusterapp/software/linux-x86_64-centos5/R-3.2.3/lib64/R/library"))

suppressMessages(library(argparse))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--regions_file", type="character",help="Path to region coordinates, format: three columns with chr, start and end. Format for chr is, for example, chr10")
parser$add_argument("--name_regions", type="character",help="Identification for the regions of interest")
parser$add_argument("--genome_library", type="character", help="Human Genome library, for example, BSgenome.Hsapiens.UCSC.hg19 or BSgenome.Hsapiens.NCBI.GRCh38")
parser$add_argument("--type_of_C", type="character", default="CpG", help="To be used in cases where the analysis is focused on nonCpG Cs, in that case --type_of_C = nonCpG")
parser$add_argument("--chr", type="character",help="if hg19 chr should be of format chr3, if GRCh38 it should be just 3, no chr")

args <- parser$parse_args()

print(args)

suppressMessages(library(args$genome_library,character.only=TRUE))

#=======================
### Auxiliary functions
#=======================

## the function CpG_coordinates_per_chr_function takes the chr and start/end coordinates of regions in a chromosome plus the genome of interest and
## find the coordinates from each CpG inside that region considering the POSITIVE STRAND
CpG_coordinates_per_chr_function <- function(region_coord,genome,type_C){

  CpGs_per_region_coordinates <- NULL

  for(i in 1:dim(region_coord)[1]){
    # print(i)

    if(type_C == "CpG"){ ### extracting coordinates of CpG Cs
      C_coords <- start(matchPattern("CG", genome[[as.character(region_coord$chr[i])]][region_coord$start[i]:region_coord$end[i]]))
    }

    if(type_C == "nonCpG"){ ### extracting coordinates of nonCpG Cs
      Cs <- as.vector(gregexpr("C", genome[[as.character(region_coord$chr[i])]][region_coord$start[i]:region_coord$end[i]])[[1]])
      CGs <- as.vector(gregexpr("CG", genome[[as.character(region_coord$chr[i])]][region_coord$start[i]:region_coord$end[i]])[[1]])
      C_coords <- setdiff(Cs,CGs)
    }

    if (length(C_coords) == 0){
      C_coords_positions <- NULL
    } else {
      C_coords_positions <- (C_coords-1) + region_coord$star[i]

    }


    if(length(C_coords_positions) == 0){
      CpGs_per_region_coordinates  <- rbind(CpGs_per_region_coordinates,NULL)
    } else {

      ### for Cecilia
      #table <- cbind(rep(region_coord$chr[i],length(CGs_positions)),CGs_positions)
      #colnames(table) <- NULL

      ### for other data should be this way:
      table <- cbind(rep(sub("chr","",region_coord$chr[i]),length(C_coords_positions)),C_coords_positions,C_coords_positions,rep(region_coord$start[i],length(C_coords_positions)),rep(region_coord$end[i],length(C_coords_positions)),
                     rep(length(C_coords_positions),length(C_coords_positions)), rep(region_coord$region_length[i],length(C_coords_positions)), rep(region_coord$region_id[i],length(C_coords_positions)))
      colnames(table) <- NULL

      CpGs_per_region_coordinates  <- rbind(CpGs_per_region_coordinates,table)

    }

  }

  return(CpGs_per_region_coordinates)

}

#========================
### END OF FUNCTIONS
#========================

if((args$genome_library == "BSgenome.Hsapiens.UCSC.hg19") | (args$genome_library == "BSgenome.Hsapiens.NCBI.GRCh38")){
  genome_used <- Hsapiens
}

if((args$genome_library == "BSgenome.Mmusculus.UCSC.mm10")){
  genome_used <- Mmusculus
}

print("Genome")
print(genome_used)
print(args$genome_library)

region_interest <- args$name_regions

region_coordinates <- read.table(args$regions_file,header=TRUE)  ### for some reason argument sep="\t" leads to a weird header --> chr.start.end

print(dim(region_coordinates))

print(head(region_coordinates))

### TO TEST
#region_coord <- region_coordinates[1:10,]

region_coord <- region_coordinates

region_length <- region_coord$end - region_coord$start

region_coord <- cbind(region_coord,region_length)

region_coord$chr <- as.character(region_coord$chr)

region_coord$chr <- as.factor(region_coord$chr)

region_id <- paste0(region_coord[,1],":",region_coord[,2],"-",region_coord[,3],sep="")

region_coord <- cbind(region_coord,region_id)
region_coord$region_id <- as.character(region_coord$region_id)

print("region coord")
print(head(region_coord))
print(dim(region_coord))

if( (args$genome_library == "BSgenome.Hsapiens.UCSC.hg19") | (args$genome_library == "BSgenome.Mmusculus.UCSC.mm10") ){

    if( sum( grepl(x=region_coord$chr, pattern="chr") ) == length(region_coord$chr) ){ ### chr column should contain chr, that is, chr1, chr2, etc
        print("right format for chr column")
    } else {
        print("putting chr column on the right format")
        region_coord$chr <- paste0("chr",region_coord$chr)
    }

}

if(args$genome_library == "BSgenome.Hsapiens.NCBI.GRCh38"){

    if( sum(grepl(x=region_coord$chr, pattern="chr")) == 0 ){ ### chr column should not contain chr, just 1,2,..., X,Y
        print("right format for chr column")
    } else {
        print("putting chr column on the right format")
        region_coord$chr <- sub(region_coord$chr,pattern="chr",replacement="")
    }

}

print(head(region_coord))
print(dim(region_coord))

if(!is.null(args$chr)){ ### it is supposed to be always by chr, the argument chr should always be non-empty

  region_coord$chr <- as.character(region_coord$chr)

  #print((sum(region_coord$chr==args$chr)))

  if(sum(region_coord$chr==args$chr) != 0 ){

    region_tmp <- region_coord[region_coord$chr==args$chr,]

    outdir <- args$output_directory

    ### output for Cecilia, all CpGs in hg19
    #cat(sapply(c("chrm", "CpGposition"), toString), file= file.path(outdir, paste0("CpGs_coordinates_",region_interest,"_",args$chr,".tsv")), sep="\t")
    #cat("\n", file=file.path(outdir, paste0("CpGs_coordinates_",region_interest,"_",args$chr,".tsv")), append=TRUE)

    cat(sapply(c("chr","CpG_start","CpG_end","region_start","region_end","region_cpgNum","region_length","region_id"), toString), file= file.path(outdir, paste0("CpGs_coordinates_",region_interest,"_",args$chr,".tsv")), sep="\t")
    cat("\n", file=file.path(outdir, paste0("CpGs_coordinates_",region_interest,"_",args$chr,".tsv")), append=TRUE)

    CpG_coordinates_per_chr  <- NULL

    CpG_coordinates_per_chr <- CpG_coordinates_per_chr_function(region_coord = region_tmp,genome=genome_used,type_C=args$type_of_C)

    if(!is.null(CpG_coordinates_per_chr)){

      write.table(CpG_coordinates_per_chr,file=file.path(outdir, paste0("CpGs_coordinates_",region_interest,"_",args$chr,".tsv")),row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE,quote=FALSE)

    }

    system(paste0("gzip --force ", file.path(outdir, paste0("CpGs_coordinates_",region_interest,"_",args$chr,".tsv"))))

  } else {
    print(paste0("Error: chromosome ",args$chr," contain no region of interest"))
  }


} else {
    print("Error: argument for chromosome needs to be specified")
}





