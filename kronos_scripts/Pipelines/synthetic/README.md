## Pipeline for generating synthetic methylation data and running Epiclomal


### What it does
* Generates synthetic methylation data
* Runs Epiclomal, several models, with many repeats
* Finds the run with the best ELBO
* Measures average accuracy (V-measure etc) over many sets with the same parameters
* Draws plots with the results

### Requirements
  * R v3.2.3
  * Python 2.7.5

### Inputs
  * sets of parameters

### Output
  * tables and figures

### Changelog
  * 

### How to Run
See documentation http://wiki.mo.bccrc.ca/display/lab/Pipeline+Runner

### Sample file format

#sample_id	num_loci	num_clones	num_cells	clone_prevalence	error_probability	missing_probability	seed	verbose	saveall	epiclomal_runs	epiclomal_model
DATA_1	1000	3	100	0.2_0.5_0.3	0.001_0.001	0.2	1	0	0	/home/mandronescu/EPI-73_test_epiclomal_synthetic/interval_file.txt	Basic-GeMM
DATA_2	2000	4	100	0.2_0.25_0.3_0.25	0.001_0.001	0.2	1	0	0	/home/mandronescu/EPI-73_test_epiclomal_synthetic/interval_file.txt	Basic-GeMM
DATA_3	3000	5	100	0.2_0.2_0.2_0.2_0.2	0.01_0.01	0.2	1	0	0	/home/mandronescu/EPI-73_test_epiclomal_synthetic/interval_file.txt	Basic-GeMM
