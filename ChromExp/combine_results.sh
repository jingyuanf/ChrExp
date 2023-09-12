#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l highp,h_data=8G,h_rt=10:00:00
#$ -pe shared 3
#$ -M fujy2038
#$ -N COMBINE_RESULTS
#$ -o /u/home/f/fujy2038/project-ernst/Project/batch/workflow/ChromExp/log
#$ -e /u/home/f/fujy2038/project-ernst/Project/batch/workflow/ChromExp/log
#$ -m bea
#$ -t 1

echo ${SGE_TASK_ID}
source /etc/profile
module load R/4.1.0-DS

############### PARAMETERS AND DIRECTORIES ################
## LOCATION OF R FILE TO RUN
combine_across_file_name="/u/home/f/fujy2038/project-ernst/Project/batch/workflow/ChromExp/combine_results.R"


## SETTING OF WORKING DIRECTORY
wd="/u/home/f/fujy2038/project-ernst/Project/"

## SETTING OF HISTONE MARK
# mark="`sed -n ${SGE_TASK_ID}p /u/home/f/fujy2038/project-ernst/Project/batch/workflow/markers/list_of_markers.txt`"
mark="H3K4me1"
# mark="H3K4me3"
## mark="H3K9ac"
# mark="H3K9me3"
# mark="H3K27ac"
# mark="H3K27me3"
# mark="H3K36me3"
## mark="H3K79me2"


## SETTING OF PROJECT NAME
# proj_name="epimap-${mark}-hg19-non-dup_test"
proj_name="epimap-${mark}-hg19-non-dup-cluster_test/epimap-${mark}-hg19-non-dup-cluster_test_BSS00007"


## SETTING OF INPUT FILENAMES
ATLAS_CHIP="/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-${mark}-hg19/FULL-ATLAS-EPIMAP-HG19-${mark}-RAW.rds"

## SETTING OF OUTPUT DIRECTORY
across_dir="./bss-predictions/across-model/"

## GETTING TRAINING AND TESTING CELLTYPES
testing_ct="/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-${mark}-hg19/all-celltypes-non-dup-${mark}.rds"

## SET SEED
seed=42

${combine_across_file_name} --wd ${wd} \
                        -c ${ATLAS_CHIP} \
                        -p ${proj_name} \
                        -r ${across_dir} \
                        --mark ${mark} \
                        --leave_one_out \
                        --chunk_size "1e+05" \
                        --num_chunks 13 \
                        --testing_ct ${testing_ct}
