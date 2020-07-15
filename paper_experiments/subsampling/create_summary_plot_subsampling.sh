
 
Rscript ../../scripts/plot_final_results_subsampling.R \
   --output_dir PLOTS_subsampling_missp_Smallwood2014 \
   --criterion DIC_LINE_ELBOW_gainthr0.05_0.02 \
   --datasets "0.7_0.8_0.9_0.95" \
   --nreplicates "5" \
   --datapath "subsampled_missp_Smallwood2014"


Rscript ../../scripts/plot_final_results_subsampling.R \
   --output_dir PLOTS_subsampling_ncells_Smallwood2014 \
   --criterion DIC_LINE_ELBOW_gainthr0.05_0.02 \
   --datasets "0.9_0.7_0.5_0.3_0.2" \
   --nreplicates "5" \
   --datapath "subsampled_ncells_Smallwood2014"

Rscript ../../scripts/plot_final_results_subsampling.R \
   --output_dir PLOTS_subsampling_ncells_Hou2016 \
   --criterion DIC_LINE_ELBOW_gainthr0.05_0.02 \
   --datasets "0.9_0.7_0.5_0.3" \
   --nreplicates "5" \
   --datapath "subsampled_ncells_Hou2016"