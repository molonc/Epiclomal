
# Using the same script as for synthetic

declare -a CRITERIA=("DIC_LINE_ELBOW_gainthr0.05_0.02")

for CRITERION in ${CRITERIA[@]}; do

  # plot the results
  VAR="subsamp_missp"
  Rscript ../../scripts/plot_final_results_synthetic.R \
    --input_dir="subsampled_missp_Smallwood2014" \
    --output_dir="${VAR}_Smallwood2014_FINAL_RESULTS_$CRITERION/" \
    --var=$VAR \
    --criterion=$CRITERION    #>> output_${VAR}.txt 2>&1
    
done


declare -a CRITERIA=("DIC_LINE_ELBOW_gainthr0.05_0.02")

for CRITERION in ${CRITERIA[@]}; do

  # plot the results
  VAR="subsamp_ncells"
  Rscript ../../scripts/plot_final_results_synthetic.R \
    --input_dir="subsampled_ncells_Smallwood2014" \
    --output_dir="${VAR}_Smallwood2014_FINAL_RESULTS_$CRITERION/" \
    --var=$VAR \
    --criterion=$CRITERION    #>> output_${VAR}.txt 2>&1
    
done

echo "RUN FINISHED!!!"

exit


 
Rscript ../../scripts/plot_final_results_subsampling.R \
   --output_dir PLOTS_subsampling_missp_Smallwood2014 \
   --criterion DIC_LINE_ELBOW_gainthr0.05_0.02 \
   --datasets "0.7_0.8_0.9_0.95" \
   --nreplicates "5" \
   --datapath "subsampled_missp_Smallwood2014"


Rscript ../../scripts/plot_final_results_subsampling.R \
   --output_dir PLOTS_subsampling_ncells_Smallwood2014 \
   --criterion DIC_LINE_ELBOW_gainthr0.05_0.02 \
   --datasets "0.9_0.8_0.7_0.6_0.5_0.4_0.3_0.2_0.1" \
   --nreplicates "5" \
   --datapath "subsampled_ncells_Smallwood2014"

Rscript ../../scripts/plot_final_results_subsampling.R \
   --output_dir PLOTS_subsampling_ncells_Hou2016 \
   --criterion DIC_LINE_ELBOW_gainthr0.05_0.02 \
   --datasets "0.9_0.7_0.5_0.3" \
   --nreplicates "5" \
   --datapath "subsampled_ncells_Hou2016"
