suppressMessages(library(argparse))
suppressMessages(library(REpiclomal))

#======================
# arguments
#======================

# create parser object
parser <- ArgumentParser()

parser$add_argument("--output_directory", type="character", help="Path to the output directory")
parser$add_argument("--methylation_file", type="character", help="Path to methylation data")
parser$add_argument("--regions_file", type="character",help="Path to region coordinates")
parser$add_argument("--use_cache", default="1", type="integer", help="If 1, use cached data if available to save on compute time, if 0, recompute data")

args <- parser$parse_args()

print(args)

outdir <- args$output_directory
dir.create(file.path(outdir), showWarnings = FALSE)
input_CpG_data_file <- args$methylation_file
input_regions_file <- args$regions_file
use_cache <- args$use_cache

#======================
# loading the data
#======================

data <- load_data(outdir, input_CpG_data_file, input_regions_file, use_cache)

input_CpG_data <- data$input_CpG_data
input_regions <- data$input_regions
mean_meth_matrix <- data$mean_meth_matrix
rm(data)

R <- dim(input_regions)[1] ## number of regions
M <- dim(input_CpG_data)[2] ## number of loci

print("number of regions:")
print(R)
print("number of loci:")
print(M)