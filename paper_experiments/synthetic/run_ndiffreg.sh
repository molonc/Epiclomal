
declare -a arr=(1 2 3 4 5 6)

# also to add 30 repeats each

for i in "${arr[@]}"
do
	# You can access them using echo "${arr[0]}", "${arr[1]}" also
    for run in {1..10}
    do
        echo "Starting ndiff reg $i, run $run"
	    snakemake -s ../../snakemake/synthetic_data/Snakefile \
  		    -p --nolock --rerun-incomplete \
      		--configfile config/config_ndiffreg.yaml \
  	    	--config "num_reg_flip=\"$i\"" "OUTPUT_DIR=\"ndiffreg/${i}_${run}\"" "DATA_DIRECTORY=\"ndiffreg/${i}_${run}\"" \
  		    --cluster 'qsub -V -hard -q shahlab.q -l h_vmem={resources.h_vmem}G -S /bin/bash -o {params.qsub_out} -e {params.qsub_err}' \
  		    -j 20 > log/snakemake_ndiffreg_${i}_${run}.log 2>&1 &
    done

done


#snakemake -s ../../snakemake/synthetic_data/Snakefile -p --configfile config/config_ndiffreg.yaml --config "num_reg_flip=\"2\"" --cluster 'qsub -V -hard -q shahlab.q -l h_vmem={resources.h_vmem}G -S /bin/bash -o {params.qsub_out} -e {params.qsub_err}' -j 100

