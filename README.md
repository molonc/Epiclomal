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
