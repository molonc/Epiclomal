# Examples on how to create synthetic data, run Epiclomal and plot the results.
This folder contains simple examples. For complex, real life examples, see kronos_scripts. 

* Run "sh generate_synthetic_data.sh" to create a synthetic data set. Edit this script with your paths and to change the input parameters. The generated files will be in synthetic/data.

* Run "sh run_non_probabilistic_methods.sh" to run the non-probabilistic methods on the synthetic data created. This will create directory simple_hclust with the results and plots. Edit this script with your paths and to change the input parameters.

* Run "sh run_epiclomal_basic.sh" to run EpiclomalBasic 50 or other number of times. Edit the script accordingly. EpiclomalBasic will use the results of the non-probabilistic methods for some of the initializations. 

* Run "sh run_eval_basic.sh" to evaluate EpiclomalBasic's results. Directory results_basic will contain the results including visualization plots for several DIC elbow parameters. Edit ../scripts/eval_epiclomal.R (bottom) to change DIC elbow parameters.

* Run "sh run_epiclomal_region.sh" to run EpiclomalRegion 50 or other number of times. Edit the script accordingly. EpiclomalRegion will use the results of the non-probabilistic methods for some of the initializations. 

* Run "sh run_eval_region.sh" to evaluate EpiclomalRegion's results. Directory results_region will contain the results including visualization plots for several DIC elbow parameters. Edit ../scripts/eval_epiclomal.R (bottom) to change DIC elbow parameters.

