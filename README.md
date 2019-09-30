# Epiclomal
Epiclomal package, software for clustering of sparse DNA methylation data

This folder contains the following:

- epiclomal is the software for clustering in Python 3.0
- examples contains some simple examples of how to generate synthetic data, run the non-probabilistic methods, run Epiclomal and evaluate and plot the results. These examples are just bash scripts and do not use the kronos pipeline.
- kronos_scripts contains pipelines and components for running kronos pipelines. You need kronos pipeliner 2.3 to run this, see https://pypi.org/project/kronos-pipeliner/.
- process_real_data contains R scripts to pre-process DNA methylation data given a set of functional regions.
- scripts contains R scripts to generate synthetic data, run non-probabilistic methods and generate plots. Some of the requirements for scripts are: MCMCpack, densityCut (https://bitbucket.org/jerry00/densitycut_dev) and its requirements, NbClust, pcaMethods, pheatmap, argparse.

If you use this software, please cite "Epiclomal: probabilistic clustering of sparse single-cell DNA methylation data,
Camila P. E. de Souza, Mirela Andronescu, Tehmina Masud, Farhia Kabeer, Justina Biele, Emma Laks, Daniel Lai, Jazmine Brimhall, Beixi Wang, Edmund Su, Tony Hui, Qi Cao, Marcus Wong, Michelle Moksa, Richard A. Moore, Martin Hirst, Samuel Aparicio, Sohrab P. Shah, doi: https://doi.org/10.1101/414482"


## Setup and Installation

Set up conda with the required packages.

First ensure you have the correct channels:
```
conda config --add channels 'https://conda.anaconda.org/dranew'
conda config --add channels 'https://conda.anaconda.org/aroth85'
conda config --add channels 'https://conda.anaconda.org/shahcompbio'
conda config --add channels 'bioconda'
conda config --add channels 'r'
conda config --add channels 'conda-forge'
```

### From Source

Clone Epiclomal:

```
git clone https://github.com/shahcompbio/Epiclomal.git
cd Epiclomal
```

Then create an environment with the required packages:

```
conda create --name Epiclomal --file conda_packages.txt
```

Activate the environment:

```
source activate Epiclomal
```

Add Epiclomal Python package into the current site packages:
```
python setup.py install
```

Add Epiclomal R package into current site packages:
```
R CMD build REpiclomal
R CMD INSTALL REpiclomal_1.0.tar.gz
```

## Usage

### Run entire pipeline with generated synthetic data

A Snakemake workflow exists to generate synthetic data and run the clustering and cluster evaluation software against the generated data.

This workflow follows this diagram, but with 300 iterations for run_epiclomal_basic and run_epiclomal_region
![Alt text](./dag.svg)

To run the Snakemake workflow, first edit the config file found at Epiclomal/snakemake/synthetic_data/config.yaml with appropriate paths and parameters. Then run
```
snakemake -s /path/to/Epiclomal/snakemake/synthetic_data/Snakefile
```
To run the workflow locally. To submit the jobs on the shahlab cluster and with parallelization, run
```
snakemake -s /path/to/Epiclomal/snakemake/synthetic_data/Snakefile --cluster 'qsub -V -hard -q shahlab.q -l h_vmem=8G -P shahlab_high -S /bin/bash' -j 32
```

To run each step of the synthetic data workflow individually, follow the steps outlined here: https://github.com/shahcompbio/Epiclomal/blob/master/examples/README.md

### Using individual components
### Epiclomal R package
In R, run `library(REpiclomal)` to use Epiclomal R Package and `?REpiclomal` for documentation.

## Input
## Output


