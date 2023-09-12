#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l highp,h_data=8G,h_rt=12:00:00
#$ -pe shared 8
#$ -M fujy2038
#$ -N Preprocess_EpiMap_rm_dup
#$ -o ./preprocess_epimap/log
#$ -e ./preprocess_epimap/log
#$ -m bea
#$ -t 1-8

echo ${SGE_TASK_ID}
source /etc/profile
module load R/4.1.0-DS

marker="`sed -n ${SGE_TASK_ID}p ./markers/list_of_markers.txt`"
run_file="./preprocess_epimap/remove_duplicates_and_quantile_normalize.R"

${run_file} -m ${marker} \
            --wig \
            --rm_dup \

