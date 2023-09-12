#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h_data=8G,h_rt=1:00:00
#$ -pe shared 4
#$ -M user_name
#$ -N Deconvolve
#$ -o ./postprocess_epimap/log
#$ -e ./postprocess_epimap/log
#$ -m bea
#$ -t 1

echo ${SGE_TASK_ID}
source /etc/profile
module load R/4.1.0-DS

postprocess_dir="./postprocess_epimap/"
run_file="${postprocess_dir}/deconvolve_pred.R"
mark="`sed -n ${SGE_TASK_ID}p /u/home/f/fujy2038/project-ernst/Project/batch/workflow/markers/list_of_markers.txt`"
chr="chr1"
proj_name="epimap-${mark}-hg19-non-dup"
dir_name_chr="epimap-${mark}-chromimpute-rna-all-pred-${chr}"
pred_dir="/u/home/f/fujy2038/project-ernst/Project/bss-predictions/"
testing_ct="/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-${mark}-hg19/all-celltypes-non-dup-${mark}.rds"
training_ct="/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-${mark}-hg19/all-celltypes-non-dup-${mark}.rds"
bulk_tracks_f="/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-${mark}-hg19/FULL-ATLAS-EPIMAP-HG19-${mark}-RAW.rds"
log_file="/u/home/f/fujy2038/project-ernst/Project/batch/workflow/postprocess_epimap/log/deconvolve_pred_${proj_name}_${chr}.Rout"

Rscript --no-save --no-restore --verbose ${run_file} ${mark} ${chr} ${proj_name} ${dir_name_chr} ${pred_dir} ${testing_ct} ${training_ct} ${bulk_tracks_f} > ${log_file} 2>&1
