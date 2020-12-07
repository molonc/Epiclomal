snakemake -s ../../../snakemake/process_real_data/Snakefile \
  --configfile config_Smallwood2014.yaml
#  --cluster 'qsub -V -hard -q shahlab.q -l h_vmem={resources.h_vmem}G -S /bin/bash -o {params.qsub_out} -e {params.qsub_err}' -j 30
