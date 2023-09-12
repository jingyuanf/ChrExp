#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l highp,h_data=8G,h_rt=2:00:00
#$ -pe shared 4
#$ -M fujy2038
#$ -N ChromImpute-Train
#$ -o ./ChromImpute/log
#$ -e ./ChromImpute/log
#$ -m bea
#$ -t 1-16

source /etc/profile
echo ${SGE_TASK_ID}
module load R/4.1.0-DS

num_cell_to_convolve=0
num_train_ct=20
num_group=0
# ct="`sed -n ${SGE_TASK_ID}p ./holdout-celltypes-rm-dup/BSS_training_validation_celltypes/testing_ct_16_celltypes_${num_group}.txt`"
### Human epimap
get_ct_list="./holdout-celltypes-rm-dup/BSS_training_validation_celltypes/training_ct_20_celltypes_0.rds"
# ct="`sed -n ${SGE_TASK_ID}p ./holdout-celltypes-rm-dup/BSS_training_validation_celltypes/validing_ct_16_celltypes_${num_train_ct}_${num_group}.txt`"
ct="`sed -n ${SGE_TASK_ID}p ./holdout-celltypes-rm-dup/BSS_training_validation_celltypes/testing_ct_16_celltypes_${num_group}.txt`"

### Mouse brain
# get_ct_list="./Data/mouse_bulk_to_sc_processed/training_celltypes_sc.rds"
# ct="`sed -n ${SGE_TASK_ID}p ./Data/mouse_bulk_to_sc_processed/testing_celltypes_sc.txt`"

Rscript --no-save --no-restore --verbose ./ChromImpute/generate_mark_tables.R ${num_cell_to_convolve} ${num_train_ct} ${num_group} ${get_ct_list} ${ct} > ./ChromImpute/log/generate_mark_tables.R_${num_train_ct}_${num_group}_${ct}.Rout 2>&1
# Rscript --no-save --no-restore --verbose ./ChromImpute/generate_mark_tables_mouse.R ${num_cell_to_convolve} ${num_train_ct} ${get_ct_list} ${ct} > ./ChromImpute/log/generate_mark_tables_mouse.R_${num_train_ct}_${ct}.Rout 2>&1