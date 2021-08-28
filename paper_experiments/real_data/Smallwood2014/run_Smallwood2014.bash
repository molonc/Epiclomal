#!/bin/bash
echo "test"
declare -a arr=("0_0.95_10000")
#				"0_1_0.01"
#				"0_1_0.05"    
#				"0_1_0.1"
#				"0_0.95_0.01"
#				"0_0.95_0.05"    
#				"0_0.95_0.1")

pwd
for i in "${arr[@]}"
do
   echo "Starting cutoffs $i"
	# You can access them using echo "${arr[0]}", "${arr[1]}" also

	snakemake -s ../../../snakemake/real_data/Snakefile \
  		-p -r --nolock \
  		--configfile config_Smallwood2014.yaml \
  		--config "CUTOFFS=\"$i\"" \
  		# --cluster 'qsub -V -hard -q shahlab.q -l h_vmem={resources.h_vmem}G -S /bin/bash -o {params.qsub_out} -e {params.qsub_err}' \
  		-j 30 > snakemake_log_$i.txt 2>&1 &

done
 
 
#snakemake -s ../../../snakemake/real_data/Snakefile --unlock --configfile config_Smallwood2014.yaml

