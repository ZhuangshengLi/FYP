# FYP
---

## Data Source

The data used in this project is from the paper **"Formative pluripotent stem cells show features of epiblast cells poised for gastrulation"**, with GEO accession number **GSE154290**. The gene dataset can be accessed at 

[https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA645812&o=acc_s%3Aa](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA645812&o=acc_s%3Aa).

Specific genome information selected can be found in the `idname` file.

---

## Required Files

### 1. HISAT2 Alignment Index File

The HISAT2 alignment index file needs to be downloaded from the HISAT2 GitHub page. The version used in this project is **mm10**, and the download link is as follows:

[https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz](https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz)

After downloading, extract the file and store it in a local directory for HISAT2 usage.

### 2. GTF Annotation File

The GTF annotation file is required for generating the gene expression count matrix using featureCounts. The version used in this project is **GRCm38**, and the download link is:

[https://www.gencodegenes.org/mouse/release_M25.html](https://www.gencodegenes.org/mouse/release_M25.html)

### 3. Salmon Index File

To build the Salmon index, the reference transcriptome sequence `cdna.fa` file must be downloaded. The file can be accessed from the Ensembl database:

[https://nov2020.archive.ensembl.org/Mus_musculus/Info/Index](https://nov2020.archive.ensembl.org/Mus_musculus/Info/Index)

After downloading, run the following commands to build the index:

```bash
mkdir ~/fyp/reference/index/salmon/grcm38
cd ~/fyp/reference/index/salmon/grcm38
salmon index -p 12 -t ~/reference/ensembl/grcm38.cdna.fa.gz -i ./
```

### 4. Sample Group Annotation File

To annotate sample groups, download the sample Metadata information table `SraRunTable.txt` from the GEO website.

![Metadata](https://github.com/ZhuangshengLi/FYP/blob/main/pic/Metadata.jpg)

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
   - Aligns sequence data using HISAT2 and generates BAM format files.

6. `4_counts_featurecounts.sh`  
   - Generates a gene expression count matrix using featureCounts.

7. `gtf_geneid2symbol_gencode.sh`  
   - Maps gene IDs to gene symbols.

### **Notes**
- **Gene Sequence Set Name**: The name of the gene sequence set must be manually input and stored in the `idname` file.
- **Gene Annotation File**: The required gene annotation file must be downloaded from the internet. The file names may differ from those in the scripts, so manual updates to the script may be needed before execution.
- **Large Files**: Due to GitHub's single file size limit (files cannot exceed 100MB), the gene annotation files are not included in this repository and must be obtained independently.

---

## R Script Functions

R scripts are located in the `test/R/` directory and serve the following purposes:

1. `Relevance.R`  
   - Generates visualizations representing the correlation between normalization methods to aid in comparing their relationships.

2. `envDep.R`  
   - Records the R packages used in data analysis to facilitate environment setup.

3. `getCountsMatric.R`  
   - Extracts count matrices from the results of different normalization methods.

4. `outputMatric.R`  
   - Generates count matrices formatted for use with various normalization methods.

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
```

---

### **Subsequent Steps**

After completing the above steps, import the `SraRunTable.txt` file and other output results into R to process the data and convert it into a gene expression count matrix for further analysis.


### Attribution
The R scripts `test/R/1_TMM_norm.r`, `test/R/1_TMM_norm_fyp.r`, and `test/R/custom_functions.r` have been adapted from the code prototypes available at ML_PDACBiomarker. These files have been modified for this project and make use of the dataset provided by the [ML_PDACBiomarker project](https://github.com/Victormah/ML_PDACBiomarker).
