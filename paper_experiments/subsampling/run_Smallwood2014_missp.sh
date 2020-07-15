
SDIR="subsampled_missp_Smallwood2014"

#declare -a arr=`ls $SDIR`



arr=($(ls -d $SDIR/*))

counter=0
for i in "${arr[@]}"
do
   echo "Starting dir $i"
	# You can access them using echo "${arr[0]}", "${arr[1]}" also

    counter=$((counter+1))
	snakemake -s ../../snakemake/real_data/Snakefile \
  		-p -r --nolock \
  		--configfile config_Smallwood2014_missp.yaml \
  		--config "CUTOFFS=\"$i\"" \
  		--cluster 'qsub -V -hard -q shahlab.q -l h_vmem={resources.h_vmem}G -S /bin/bash -o {params.qsub_out} -e {params.qsub_err}' \
  		-j 500 > snakemake_missp_$counter.txt 2>&1
    #echo "Done"
    #exit
done
 
 
#snakemake -s ../../../snakemake/real_data/Snakefile --unlock --configfile config_Smallwood2014.yaml

