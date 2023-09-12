#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l highp,h_data=8G,h_rt=24:00:00
#$ -pe shared 8
#$ -M fujy2038
#$ -N Preprocess_EpiMap_gene_exp
#$ -o ./preprocess_epimap/log
#$ -e ./preprocess_epimap/log
#$ -m bea

# echo ${SGE_TASK_ID}
source /etc/profile
module load R/4.1.0-DS

# marker="`sed -n ${SGE_TASK_ID}p ./markers/list_of_markers.txt`"
run_file="./preprocess_epimap/preprocess_atlas_gene.R"
wd="./mypath"
gene_train="./data/epimap-gene-exp/FULL-ATLAS-EPIMAP-HG19-RNA-GENEID-LOG2FPKM.rds"
gene_test="./data/epimap-gene-exp/FULL-ATLAS-EPIMAP-HG19-RNA-GENEID-LOG2FPKM.rds"
gene_meta="./data/epimap-metadata/gencode_gene_metadata.rds"
new_file_location="epimap-gene-exp"

Rscript --no-save --no-restore --verbose ${run_file} ${wd} ${gene_train} ${gene_test} ${gene_meta} ${new_file_location}

