#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l highp,h_data=4G,h_rt=3:00:00
#$ -pe shared 6
#$ -M fujy2038
#$ -N ChromImpute-Train
#$ -o ./ChromImpute/log
#$ -e ./ChromImpute/log
#$ -m n
#$ -t 12

source /etc/profile
echo ${SGE_TASK_ID}


### 8 histone marks + 1 RNA
# mark="`sed -n ${SGE_TASK_ID}p ./markers/list_of_markers_w_rna.txt`"
mark_rna="RNA"
mark="H3K4me1"

### Full mark table (used in first step for converting tracks)
# inputinfofile="./bss-predictions/chromimpute/epimap/mark_table_full_gb.txt"

### Parallelize over sample
# sample="`sed -n ${SGE_TASK_ID}p ./Correct_data/epimap-gene-exp/gene_exp_celltypes.txt`"
sample="`sed -n ${SGE_TASK_ID}p ./Correct_data/epimap-${mark}-hg19/all-celltypes-non-dup-${mark}.txt`"

### Individual mark table (separate table for each mark, used in following steps)
inputinfofile="./bss-predictions/chromimpute/epimap/mark_tables/mark_table_${mark}.txt"
chrominfofile="./bss-predictions/chromimpute/epimap/chrom_info_chr1.txt"
# chrominfofile="./bss-predictions/chromimpute/epimap/chrom_info.txt"

### Mark table for RNA
# inputinfofile="./bss-predictions/chromimpute/epimap/mark_tables/mark_table_${mark_rna}.txt"


### Input directory of histone modification tracks
# INPUTDIR="./Correct_data/epimap-${mark}-hg19/raw-tracks-wig"

### Input directory of RNA tracks
INPUTDIR="./Correct_data/epimap-gene-exp/raw-tracks-genebody/"


CONVERTEDDIR="./bss-predictions/chromimpute/epimap/converted_tracks/"
DISTANCEDIR="./bss-predictions/chromimpute/epimap-${mark}-chromimpute-rna-all-pred-chr1/distance/"
TRAINDATADIR="./bss-predictions/chromimpute/epimap-${mark}-chromimpute-rna-all-pred-chr1/train_data/"
PREDICTORDIR="./bss-predictions/chromimpute/epimap-${mark}-chromimpute-rna-all-pred-chr1/predictors/"
OUTPUTIMPUTEDIR="./bss-predictions/chromimpute/epimap-${mark}-chromimpute-rna-all-pred-chr1/output/"


mkdir -p ${PREDICTORDIR}
mkdir -p ${OUTPUTIMPUTEDIR}

java -jar ./ChromImpute/ChromImpute.jar Convert -r 200 -l ${sample} $INPUTDIR $inputinfofile $chrominfofile $CONVERTEDDIR
java -jar ./ChromImpute/ChromImpute.jar ComputeGlobalDist -r 200 -m ${mark_rna} $CONVERTEDDIR $inputinfofile $chrominfofile $DISTANCEDIR
java -jar ./ChromImpute/ChromImpute.jar GenerateTrainData -r 200 -tieglobal $CONVERTEDDIR $DISTANCEDIR $inputinfofile $chrominfofile $TRAINDATADIR $mark
java -jar ./ChromImpute/ChromImpute.jar Train $TRAINDATADIR $inputinfofile $PREDICTORDIR $sample $mark
java -jar ./ChromImpute/ChromImpute.jar Apply -r 200 -tieglobal $CONVERTEDDIR $DISTANCEDIR $PREDICTORDIR $inputinfofile $chrominfofile $OUTPUTIMPUTEDIR $sample $mark


