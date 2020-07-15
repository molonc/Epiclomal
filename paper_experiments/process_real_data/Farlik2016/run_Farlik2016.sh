snakemake -s ../../../snakemake/process_real_data/Snakefile \
  -p \
  --configfile config_Farlik2016.yaml \
  --cluster 'qsub -V -hard -q shahlab.q -l h_vmem={resources.h_vmem}G -S /bin/bash -o {params.qsub_out} -e {params.qsub_err}' \
  -j 30
