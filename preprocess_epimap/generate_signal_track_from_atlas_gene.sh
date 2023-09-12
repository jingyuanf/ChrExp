#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l highp,h_data=8G,h_rt=1:00:00
#$ -pe shared 4
#$ -M user_name
#$ -N Preprocess_EpiMap_generate_exp_track
#$ -o ./preprocess_epimap/log
#$ -e ./preprocess_epimap/log
#$ -m bea
#$ -t 1-340

# echo ${SGE_TASK_ID}
# Total 340 cell types
source /etc/profile
module load R/4.1.0-DS

# marker="`sed -n ${SGE_TASK_ID}p ./markers/list_of_markers.txt`"
run_file="./preprocess_epimap/generate_signal_track_from_atlas_gene.R"
ct="`sed -n ${SGE_TASK_ID}p ./data/epimap-gene-exp/gene_exp_celltypes.txt`"

${run_file} --ct ${ct} \

