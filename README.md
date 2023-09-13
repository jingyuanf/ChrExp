ChromExp is a regression-based ensemble method that predicts chromatin signal from gene expression. This is the pipeline for evaluating ChromExp model and ChromImpute performance on EpiMap data. Notice that the example job submission code (.sh) files are based on the UCLA Hoffman2 Linux compute cluster. For other platforms, please adjust the submission code accordingly.

## Step 0: Specify your project directory

Specify your project directory here (my example project directory is listed here):
`proj_dir="/u/project/fujy2038/Project"` 

Note that for some R files, the filenames and project directory are currently hardcoded into the file. You might need to change them while applying them to other files.

## Step 1: Preprocessing bulk gene expression and chromatin signal data from EpiMap


### 1.1 Download Gene Expression data from EpiMap

Website: https://personal.broadinstitute.org/cboix/epimap/rnaseq_data/

Download "merged_log2fpkm.mtx.gz" from website (there are also quantile normalized version, "merged_qn_log2fpkm.mtx.gz"):

Set up directory for storing gene expression data:
`mkdir -p ${proj_dir}/Correct_data/epimap-gene-exp/`

Download RNA-seq from EpiMap website:
`wget -P ${proj_dir}/Correct_data/epimap-gene-exp/ https://personal.broadinstitute.org/cboix/epimap/rnaseq_data/merged_log2fpkm.mtx.gz`

Quantile Normalized version:
`wget -P ${proj_dir}/Correct_data/epimap-gene-exp/ https://personal.broadinstitute.org/cboix/epimap/rnaseq_data/merged_qn_log2fpkm.mtx.gz`

### 1.2 Download Gene Annotation Metadata from ENCODE

Website: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

Download "gencode.v19.annotation.gtf.gz" from website:

Set up directory for storing gene annotation metadata:
`mkdir -p ${proj_dir}/Correct_data/epimap-metadata/`

Download Gene Annotation Metadata from ENCODE:
`wget -P ${proj_dir}/Correct_data/epimap-metadata/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz`

### 1.3 Process Gene Expression and Gene Annotation data

1.3.1 Run /preprocess_epimap/generate_metadata.R to process gene annotation meta data
Meta data are subset to genes, protein-coding genes, chr1-X, KNOWN genes

1.3.2 Run /preprocess_epimap/generate_atlas_gene.R to process gene expression data

### 1.4 Download Chromatin Signal data from EpiMap

We download the chromatin signal data from EpiMap from this website: 
https://epigenome.wustl.edu/epimap/data/observed/

Run "preprocess_epimap/download_from_epimap_and_merge.sh". Right now it is processing the data by making 200 bp bins. 

### 1.5 Remove duplicated cell types
EpiMap Metadata is here: https://personal.broadinstitute.org/cboix/epimap/metadata/Imputation_Metadata.xlsx

I've downloaded it to "./data/Imputation_Metadata.xlsx"

Run "/preprocess_epimap/remove_duplicates_and_quantile_normalize.sh"

### 1.6 Standardize the format of Gene Expression and Gene Annotation data for putting into the model

Run "preprocess_epimap/preprocess_atlas_gene.sh". Right now it is processing the data by making 200 bp bins. 

## Step 2: Run ChrExp model (Gene expression -> Chromatin Signal)

Run "ChromExp/run_models.sh", run "train_model.R" first and then run "test_model.R" next.

Run "ChromExp/combine_results.sh", to combine the results from multiple chunks and from different cell types.

## Step 3: Run ChromImpute model (Gene expression -> Chromatin Signal)

### 3.1 Prepare histone mark tracks for running ChromImpute model

Run "preprocess_epimap/remove_duplicates_and_quantile_normalize.sh" with "--wig" handle to produce wig outputs.

### 3.2 Prepare gene expression tracks for running ChromImpute model

Run "preprocess_epimap/generate_signal_track_from_gene_exp.sh" to generate tracks from gene expression data. All cell types with gene expression available are located at "./data/gene_exp_celltypes.txt"

Two options: 
1. Use "--tss" option to annotate the tracks with only TSS sites having gene expression values
2. Not use "--tss" option so that all gene body are annotated with gene expression values

### 3.3 Prepare mark tables and sequence info table

1. Use "ChromImpute/generate_seqlength_table.R" to generate a sequence info table. (Sequence info file is located at "./data/chrom_info.txt")
2. Use "ChromImpute/generate_full_mark_table.sh" to generate a full mark table.
3. Then run "ChromImpute/generate_individual_mark_tables.sh" to produce mark tables specifically for each imputation task.

## Step 4: Run ChromImpute model (Chromatin Signal -> Gene expression)

Run "ChromImpute/batch_chromimpute.sh"
TODO: Note that for "ComputeGlobalDist", I only computed distance based on RNA marker (using -m command). Output distance files stored separately for each mark.

TODO: Note that "GenerateTrainData" currently uses "-c chr1" which generates train data for chromosome 1 only. But the predictions are made over all chromosomes. Later will need to make predictions on whole genome using whole genome train data

Run "ChromImpute/wigtobw.sh" to convert wig outputs to bw outputs.

## Step 5: Post-process

Run "postprocess_epimap/deconvolve_pred.sh" to generate formatted outputs and deconvolved outputs.



