#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h_data=4G,h_rt=5:00:00
#$ -pe shared 6
#$ -M fujy2038
#$ -N ChromImpute_wigtobw
#$ -o ./ChromImpute/log
#$ -e ./ChromImpute/log
#$ -m bea
#$ -t 1

echo ${SGE_TASK_ID}
source /etc/profile
module load R/4.1.0-DS

run_file="./ChromImpute/wigtobw.R"
mark="`sed -n ${SGE_TASK_ID}p ./markers/list_of_markers.txt`"
proj_name="epimap-${mark}-chromimpute-rna-all-pred-chr1"
chr="chr1"
testing_ct="./data/epimap-${mark}-hg19/all-celltypes-non-dup-${mark}.rds"
chrimp_dir="./bss-predictions/chromimpute"
log_file="./ChromImpute/log/wigtobw_${mark}_${chr}.Rout"

Rscript --no-save --no-restore --verbose ${run_file} ${mark} ${proj_name} ${chr} ${testing_ct} ${chrimp_dir} > ${log_file} 2>&1

