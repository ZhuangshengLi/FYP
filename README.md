# FYP
---

## Shell Script Execution Order

The scripts are located in the `test/sh/` directory and should be executed in the following order:

1. `00_prefetch.sh`  
   - Downloads SRA data.

2. `01_sra2fq_qc1.sh`  
   - Converts SRA data to FASTQ format and performs initial quality control.

3. `2_cleanfq_qc2.sh`  
   - Cleans FASTQ data and performs secondary quality control.

4. `3333_salmon.sh`  
   - Uses Salmon for rapid transcriptome quantification.

5. `3_align2sam2bam_hisat2.sh`  
   - Aligns sequence data using HISAT2 and generates BAM format.

6. `4_counts_featurecounts.sh`  
   - Generates a gene expression counts matrix using featureCounts.

7. `gtf_geneid2symbol_gencode.sh`  
   - Maps gene IDs to gene symbols.

### **Notes**
- **Gene Sequence Set Name**: The name of the gene sequence set must be manually input and stored in the `idname` file.
- **Gene Annotation File**: The required gene annotation file must be downloaded from the internet. File names may differ from those in the scripts, so you need to manually update the script before execution.
- **Large Files**: Due to GitHub's file size limit (single files cannot exceed 100MB), gene annotation files are not included in this repository and must be obtained independently.

---

## R Script Functions

R scripts are located in the `test/R/` directory and serve the following purposes:

1. `Relevance.R`  
   - Generates visualizations representing the correlation between normalization methods to aid in comparing their relationships.

2. `envDep.R`  
   - Records the R packages used in data analysis to facilitate environment setup.

3. `getCountsMatric.R`  
   - Extracts counts matrices from the results of different normalization methods.

4. `outputMatric.R`  
   - Generates counts matrices formatted for use with various normalization methods.

---

## Usage Instructions

### **Running Shell Scripts**
Execute the scripts in the `test/sh/` directory sequentially using the following commands:
```bash
bash 00_prefetch.sh
bash 01_sra2fq_qc1.sh
bash 2_cleanfq_qc2.sh
bash 3333_salmon.sh
bash 3_align2sam2bam_hisat2.sh
bash 4_counts_featurecounts.sh
bash gtf_geneid2symbol_gencode.sh

