#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l highp,h_data=8G,h_rt=24:00:00
#$ -pe shared 8
#$ -M user_name
#$ -N Preprocess_EpiMap_download
#$ -o ./preprocess_epimap/log
#$ -e ./preprocess_epimap/log
#$ -m bea
#$ -t 1

echo ${SGE_TASK_ID}
source /etc/profile
module load R/4.1.0-DS

marker="`sed -n ${SGE_TASK_ID}p ./markers/list_of_markers.txt`"
run_file="./preprocess_epimap/download_from_epimap_and_merge.R"

${run_file} -m ${marker}

