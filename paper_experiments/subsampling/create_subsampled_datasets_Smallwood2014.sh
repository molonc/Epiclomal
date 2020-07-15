Rscript  \
    ../../process_real_data/subsampling_ncells.R \
    --output_directory subsampled_ncells_Hou2016 \
    --data_ID Hou2016 \
    --file_input_Epiclomal ../process_real_data/Hou2016/process_Hou2016/epiclomal_input/0_0.95_10000/input_Epiclomal_process_Hou2016.tsv.gz \
    --file_final_regions ../process_real_data/Hou2016/process_Hou2016/epiclomal_input/0_0.95_10000/regionIDs_input_Epiclomal_process_Hou2016.tsv.gz \
    --file_true_membership ../process_real_data/Hou2016/data_Hou2016/true_clone_membership.txt.gz \
    --keep_cells_prop 0.9_0.7_0.5_0.3 \
    --number_subsamples 5
    
exit    

Rscript  \
    ../../process_real_data/subsampling_ncells.R \
    --output_directory subsampled_ncells_Smallwood2014 \
    --data_ID Smallwood2014 \
    --file_input_Epiclomal ../process_real_data/Smallwood2014/process_Smallwood2014/epiclomal_input/0_0.95_10000/input_Epiclomal_process_Smallwood2014.tsv.gz \
    --file_final_regions ../process_real_data/Smallwood2014/process_Smallwood2014/epiclomal_input/0_0.95_10000/regionIDs_input_Epiclomal_process_Smallwood2014.tsv.gz \
    --file_true_membership ../process_real_data/Smallwood2014/data_Smallwood2014/true_clone_membership.txt.gz \
    --keep_cells_prop 0.9_0.7_0.5_0.3_0.2_0.1 \
    --number_subsamples 5

exit

Rscript  \
    ../../process_real_data/subsampling_missp.R \
    --output_directory subsampled_missp_Smallwood2014 \
    --data_ID Smallwood2014 \
    --file_input_Epiclomal ../process_real_data/Smallwood2014/process_Smallwood2014/epiclomal_input/0_0.95_10000/input_Epiclomal_process_Smallwood2014.tsv.gz \
    --file_final_regions ../process_real_data/Smallwood2014/process_Smallwood2014/epiclomal_input/0_0.95_10000/regionIDs_input_Epiclomal_process_Smallwood2014.tsv.gz \
    --new_miss_prop 0.7_0.8_0.9_0.95 \
    --number_subsamples 5

