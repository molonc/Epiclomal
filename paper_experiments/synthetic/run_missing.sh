snakemake \
    -s ../../snakemake/synthetic_data/Snakefile \
    --cluster 'qsub -V -hard -q shahlab.q -l h_vmem={resources.h_vmem}G -S /bin/bash -o {params.qsub_out} -e {params.qsub_err}' \
    -j 8 \
    --configfile config/config_missing.yaml

