# Scripts for pre-processing the real data
These are the scripts used by the pre-process pipeline, given a set of region coordinates. 

## Task 1: CpG_coordinates_accross_regions
Extracts the CpG coordinates for the given regions, outputs one file per chromosome.
Output file example:
```
chr     CpG_start       CpG_end region_start    region_end      region_cpgNum   region_length   region_id
1       3531625 3531625 3531624 3531843 27      219     chr1:3531624-3531843
1       3531627 3531627 3531624 3531843 27      219     chr1:3531624-3531843
1       3531633 3531633 3531624 3531843 27      219     chr1:3531624-3531843
1       3531651 3531651 3531624 3531843 27      219     chr1:3531624-3531843
```

## Task 2: cell_based_methylation_extraction.R
Extracts the CpG coordinates in the given cells, outputs one file per cell (CpG level).
This output file contains all the data, before filtering, in the given set of regions, including missing data (NA) and non-binary methylation calls (column meth_frac). Output file example:
```
chr     CpG_start       CpG_end region_start    region_end      region_cpgNum   region_length   region_id       meth_frac       count_meth
      count_unmeth    cell_id
1       3531625 3531626 3531624 3531843 27      219     chr1:3531624-3531843    NA      NA      NA      GSM1370535_2i_1
1       3531627 3531628 3531624 3531843 27      219     chr1:3531624-3531843    NA      NA      NA      GSM1370535_2i_1
1       3531633 3531634 3531624 3531843 27      219     chr1:3531624-3531843    0       0       1       GSM1370535_2i_1
1       3531651 3531652 3531624 3531843 27      219     chr1:3531624-3531843    1       1       0       GSM1370535_2i_1
```

## Task 3: cell_based_stats_methylation.R
Calculates various statistics for each region in a cell, outputs one file per cell (region level). Output file example:
```
region_id       region_mean     region_median   region_IQR      region_miss     region_cov      cell_id
chr1:6214430-6215332    0       0       0       0.644859813084112       1.47368421052632        GSM1370535_2i_1
chr1:6359082-6359442    0       0       0       0.0689655172413793      1.88888888888889        GSM1370535_2i_1
chr1:6382891-6383203    0       0       0       0.238095238095238       3.125   GSM1370535_2i_1
chr1:6441200-6441612    0       0       0       0.661290322580645       2.23809523809524        GSM1370535_2i_1
```

## Task 4: stats_methylation.R
Calculates overall stats in all cells and final regions for Epiclomal

## Task 5: filter_regions.R
Finds and removes redundant regions to reduce noise of epiclomal input files. This is work in progress and was not included in the paper.

## Task 6: get_data_ready_epiclomal.R 
Creates methylation and regions files for Epiclomal's input

## Other files in folder:
- filter_regions.cpp - function to find redundant regions in filtering step is in C++ for performance optimization.