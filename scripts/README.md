# Various scripts for running non-probabilistic models, plotting and visualization

This folder contains the following:

- densitycut.R is the main script in R that runs the non-probabilistic model DensityCut.
- euclidenclust.R is the main script in R that runs the non-probabilistic model EuclideanClust.
- hammingclust.R is the main script in R that runs the non-probabilistic model HammingClust.
- pearsonclust.R is the main script in R that runs the non-probabilistic model PearsonClust.
- hclust.R implements in R the non-probabilistic models EuclideanClust, PearsonClust and HammingClust and DensityCut sequentially, then joins the results.
- generate_synthetic.R is the main script in R that creates a synthetic data set with the given parameters.
- plot_final_results.R creates final plots for the synthetic data.
- plot_final_results_real.R creates final plots for the real data.
- plot_functions.R are used by the two scripts above to plot.
