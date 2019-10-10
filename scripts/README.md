# Various scripts for running non-probabilistic models, plotting and visualization

This folder contains the following:

- ch_index.R - script called by hclust.R to apply the automatic CH index of finding the optimal number of clusters.
- hclust.R implements in R the non-probabilistic models EuclideanClust, PearsonClust and HammingClust and DensityCut. The script takes in an argument `--method` which specifies which method to run. To load and cache data, set method to `None`. Set method to `all` to run all methods sequentially, then join the results.
- generate_synthetic.R is the main script in R that creates a synthetic data set with the given parameters.
- plot_final_results.R creates final plots for the synthetic data.
- plot_final_results_real.R creates final plots for the real data.
- plot_functions.R are used by the two scripts above to plot.
- eval_epiclomal.R finds the best clustering from Epiclomal results and plots results.