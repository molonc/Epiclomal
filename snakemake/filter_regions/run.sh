snakemake -s /shahlab/pye/projects/epiclomal/Epiclomal/snakemake/filter_regions/Snakefile --configfile /shahlab/pye/projects/epiclomal/Epiclomal/snakemake/filter_regions/config.yaml --cluster 'qsub -V -hard -q shahlab.q -l h_vmem={resources.h_vmem}G -S /bin/bash' -j 32
