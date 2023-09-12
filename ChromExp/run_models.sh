#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l highp,h_data=8G,h_rt=10:00:00,highp
#$ -pe shared 4
#$ -M user_name
#$ -N MODELTRAIN
#$ -o ./ChromExp/log
#$ -e ./ChromExp/log
#$ -m bea
#$ -t 1-112

echo ${SGE_TASK_ID}
source /etc/profile
module load R/4.1.0-DS

############### PARAMETERS AND DIRECTORIES ################
## LOCATION OF R FILE TO RUN
code_dir="./ChromExp/"
train_across_file_name="${code_dir}/train_across_model_w_gene_exp_sampled_reg_tree.R"
test_across_file_name="${code_dir}/test_across_model_w_gene_exp_sampled_reg_tree.R"
clean_intermediate_file_name="${code_dir}/clean_intermediate.R"
combine_across_file_name="${code_dir}/combine_results.R"

## SETTING OF EXPERIMENT
# chr="chr${SGE_TASK_ID}"
train_chr="chr1"
test_chr="chr1"

## SETTING OF WORKING DIRECTORY
wd="./working_dir/"

## SETTING OF HISTONE MARK
mark="H3K4me1"
# mark="H3K4me3"
# mark="H3K9ac"
# mark="H3K9me3"
# mark="H3K27ac"
# mark="H3K27me3"
# mark="H3K36me3"
# mark="H3K79me2"

## SETTING OF INPUT FILENAMES
data_path="./data/"
ATLAS_GENE_train="${data_path}/epimap-gene-exp/ATLAS_GENE_train.rds"
ATLAS_GENE_test="${data_path}/epimap-gene-exp/ATLAS_GENE_test.rds"
ATLAS_CHIP="${data_path}/epimap-${mark}-hg19/raw-tracks/"
ATLAS_GENE_METADATA="${data_path}/epimap-gene-exp/ATLAS_GENE_METADATA.rds"


## SETTING OF PROJECT NAME
proj_name="epimap-${mark}-hg19"

## SETTING OF OUTPUT DIRECTORY
across_dir="./bss-predictions/across-model/"
within_dir="./bss-predictions/within-model/"
integrate_dir="./bss-predictions/integration/"

## GETTING TRAINING AND TESTING CELLTYPES
training_ct="${data_path}/epimap-${mark}-hg19/all-celltypes-non-dup-${mark}.rds"
testing_ct="`sed -n ${SGE_TASK_ID}p ${data_path}/epimap-${mark}-hg19/all-celltypes-non-dup-${mark}.txt`"

## SET SEED
seed=42

################# RUN MODELS ###################################
# ACROSS CELL-TYPE MODEL
# SLOW, RUN DIFFERENT CHROMOSOMES IN PARALLEL
dist_method="corr" ## Another option is "euc" which uses Euclidean distance to calculate nearest cell types
proj_name_train="${proj_name}_train/${proj_name}_train_${testing_ct}"
proj_name_test="${proj_name}_test/${proj_name}_test_${testing_ct}"

${train_across_file_name} -v \
                        --chr ${train_chr} \
                        --num_nearest_celltypes 5 \
                        --num_nearest_gene 5 \
			--training_ct ${training_ct} \
                        --holdout ${testing_ct} \
                        --wd ${wd} \
                        --gene_train ${ATLAS_GENE_train} \
                        -c ${ATLAS_CHIP} \
                        -m ${ATLAS_GENE_METADATA} \
                        --proj_name_train ${proj_name_train} \
                        -r ${across_dir} \
                        --seed ${seed} \
                        --dist ${dist_method} \
                        --sample_size 1000 \
                        --task_id ${SGE_TASK_ID} \
                        --mark ${mark} \
                        
                        # -k 100 \
                

for chunk in 1 2 3 4 5 6 7 8 9 10 11 12 13
do                     
        ${test_across_file_name} -v \
                                --train_chr ${train_chr} \
                                --test_chr ${test_chr} \
                                --num_nearest_celltypes 10 \
                                --num_nearest_gene 5 \
                                --training_ct ${training_ct} \
                                --testing_ct ${testing_ct} \
                                --wd ${wd} \
                                --gene_train ${ATLAS_GENE_train} \
                                --gene_test ${ATLAS_GENE_test} \
                                -c ${ATLAS_CHIP} \
                                -m ${ATLAS_GENE_METADATA} \
                                --proj_name_train ${proj_name_train} \
                                --proj_name_test ${proj_name_test} \
                                -r ${across_dir} \
                                --seed ${seed} \
                                --dist ${dist_method} \
                                --which_chunk ${chunk} \
                                --chunk_size 100000 \
                                --sample_size "1e+05" \
                                --task_id ${SGE_TASK_ID} \
                                --mark ${mark} \
                                --leave_one_out \
                                --save_intermediate \
                                -k 100 \

done


${clean_intermediate_file_name} --wd ${wd} \
                        --proj_name_test ${proj_name_test} \
                        -r ${across_dir} \

${combine_across_file_name} --wd ${wd} \
                        -c ${ATLAS_CHIP} \
                        -m ${ATLAS_GENE_METADATA} \
                        -p ${proj_name} \
                        -r ${across_dir} \

