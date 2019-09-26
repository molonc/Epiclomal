
#======================
# libraries
#======================

suppressMessages(library(argparse))

# Renaming the simple methods
# hclust -> EuclideanClust
# PBALclust -> HammingClust
# Pearsonclust -> Pearsonclust
# densitycut -> DensityCut

suppressMessages(library(REpiclomal))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--methylation_file", type="character", help="Path to methylation data")
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--max_k", type="integer",default=5, help="maximum number of clusters to be considered when cutting the tree")

parser$add_argument("--index", type="character",default="ch",help="Index to be used to choose number of clusters, default = ch")

# MA: 30Jan 2019: adding an optional imputation step
parser$add_argument("--impute", default="0", type="integer",help="If it is 1, impute with the average per region/locus, if it is 0 do nothing.")

parser$add_argument("--use_cache", default="1", type="integer", help="If 1, use cached data if available to save on compute time, if 0, recompute data")

args <- parser$parse_args()

print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)
input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file
index_type <- args$index
impute <- args$impute
use_cache <- args$use_cache

#======================
# loading the data
#======================

# Methylation data
data <- load_data(input_CpG_data_file, input_regions_file, use_cache)

input_CpG_data <- data$input_CpG_data
input_regions <- data$input_regions
rm(data)

R <- dim(input_regions)[1] ## number of regions
M <- dim(input_CpG_data)[2] ## number of loci

print("number of regions:")
print(R)
print("number of loci:")
print(M)

Max_K <- min((dim(input_CpG_data)[1]-1),args$max_k)
print(Max_K)

#print(R)
#print(M)

###################################################
### HammingClust: Tony's (PBAL manuscript) clustering approach ##
###################################################
possible_clusters <- hamming.clust(input_CpG_data, Max_K, index_type, impute, use_cache, outdir)

print("Done Hamming Clust!")




