# Scripts for pre-processing the real data
These are the scripts used by the pre-process pipeline.

- Task 1: CpG_coordinates_accross_regions - Extracts the CpG coordinates for the given regions, outputs one file per chromosome.
- Task 2: cell_based_methylation_extraction.R - Extracts the CpG coordinates in the given cells, outputs one file per cell.
- Task 3: cell_based_stats_methylation.R - Calculates various statistics for each region in a cell, outputs one file per cell.
- Task 4: stats_methylation.R - Calculates overall stats in all cells and final regions for Epiclomal
- Task 5: get_data_ready_epiclomal.R - Creates methylation and regions files for Epiclomal's input

Other scripts in this directory:
- finding_union_of_overlapping_regions.R finds the union of overlapping regions
- subsampling.R can be used to downsample a given data set for missing probability etc. This was not used in the manuscript.
