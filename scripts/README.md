# Various scripts for running non-probabilistic models, plotting and visualization

This folder contains the following:

- ch_index.R - script called by hclust.R to apply the automatic CH index of finding the optimal number of clusters.
- densitycut.R is the main script in R that runs the non-probabilistic model DensityCut.
- generate_synthetic.R is the main script in R that creates a synthetic data set with the given parameters.
- hamming_distance.R calculates the hamming distance measure. This script also does the adjustment of the Epiclomal inferred states, and calculates the naive imputation.
- hclust.R implements in R the non-probabilistic models EuclideanClust, PearsonClust and HammingClust and also calls densitycut.R.
- plot_final_results.R creates final plots for the synthetic data.
- plot_final_results_real.R creates final plots for the real data.
- plot_functions.R are used by the two scripts above to plot.
- visualization.R creates heatmaps of the methylation data and clustering result.


