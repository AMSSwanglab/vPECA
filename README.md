# vPECA 
vPECA is a Variants interpretation method by Paired Expression and Chromatin Accessibility data which can identify active and active selected regulatory elements and their gene regulatory network.


# 1. Data availability

There are three essential components for the data organization in "Chromatin accessibility landscape and regulatory network of high-altitude hypoxia adaptation".

# (1) Metadata spreadsheet
The Metadata spreadsheet is provide in Summay_info_of_HUVEC_data.xlsx, which contains four subsheets, i.e., NGS, RNA-seq, ATAC-seq, and Hi-C. Metadata refers to descriptive information about the overall study, individual samples, all protocols, and references to processed and raw data file names. Information is supplied by completing all fields of a metadata template spreadsheet. Guidelines on the content of each field are provided within the spreadsheet. The reference assembly used (hg19) is provided in the metadata spreadsheet.

# (2) Raw data
The data generated in this study, including RNA-seq, ATAC-seq, and Hi–C data were deposited at [http://www.ncbi.nlm.nih.gov/geo/] (accession number: GSE145774) and Genome Sequence Achieve (project number: CRA002025; [https://bigd.big.ac.cn/gsa/browse/CRA002025]).

The metadata, peak, and bigwig track files are available at
ftp://download.big.ac.cn/gsa/CRA002025/processed_data/.
Metadata spreadsheets provide friendly access by linking the sample names with raw and processed data download links.


# (3) Processed data
The final processed data are defined as the data on which the conclusions in the related manuscript are based. 
1) The normalized abundance measurements output from Cufflinks for the expression profiling analysis results on quantitative data for genes. 
2) Peak files with quantitative openness data with a format bed and txt files for ATAC-Seq data. 
3) Tag density files for ATAC-seq in bigWig format. They can be visualized in the UCSC genome browser. 
4) Features (e.g., genes, transcripts) in processed data files are traceable using public accession numbers or chromosome coordinates. 
5) A description of the format and content of processed data files are provided in the metadata spreadsheet data processing fields. 
6) HiC matrix file.

# (4) Data share of differetial analysis

1) Differential expressed genes
2) Differential open regions

between wildtype and adapted population and adjacent time points. Please check the detailed Note.txt file.


# 2. vPECA Source code

vPECA: variants interpretation method by Paired Expression and Chromatin Accessibility data

Version 1.0
Last updated: Jan 18, 2020

# Reference

Jingxue Xin, Hui Zhang, Yaoxi He, Zhana Duren, Caijuan Bai, Lang Chen, Xin Luo, Dong-Sheng Yan, Chaoyu Zhang, Xiang Zhu, Qiuyue Yuan, Zhanying Feng, Chaoying Cui, Xuebin Qi, Ouzhuluobu, Wing Hung Wong, Yong Wang & Bing Su. Chromatin accessibility landscape and regulatory network of high-altitude hypoxia adaptation. Nature communications 11.1 (2020): 1-20. DOI: https://doi.org/10.1038/s41467-020-18638-8. 

# Method

We develop a new method called vPECA (Variants interpretation model by Paired Expression and Chromatin Accessibility data) to model genome-wide chromatin accessibility profiles for high-altitude hypoxia adaptation in HUVEC, to reveal causal SNPs, active and active selected regulatory elements for a certain gene. vPECA can integrate our measured paired expression and chromatin accessibility data with the available public data, including population genetics data, functional genomics data in ENCODE, and Hi-C data for HUVEC. Our previous work PECA integrates paired expression and chromatin accessibility data across diverse cellular contexts and model the localization to REs of chromatin regulators (CR), the activation of REs due to CRs that are localized to them, and the effect of TFs bound to activated REs on the transcription of target genes (TG) 18. Our innovation here is to extend PECA to interpret genetic variants from population genetics and matched WGS data. vPECA models how positively selected noncoding SNPs affects the RE’s selection status, chromatin accessibility, and activity and further determine the target gene expression. The statistical modelling allows us to systematically identify active REs, active selected REs, and gene regulatory network to interpret variants.

# Processing data
vPECA model requires input as sample matched time-series RNA-seq, and ATAC-seq, and individual matched DNA-seq data together with selection scores calculated for each SNP from public data. For RNA-seq and ATAC-seq data, first we processed raw reads into an expression matrix with row genes and column samples. And chromatin accessibility data as a matrix with element by sample dimensions. The candidate RE and TG pairs based on distance are collected into a Element_gene in Data_prior.mat. Then the SNPs locate on REs and their corresponding selection scores are in a text file named element_SNP_use.txt. The prior (TF-TG, TG-RE) learned from public data could be set to certain number if it is not available. TF binding strength are calculated from motif scan algorithm.

# Running vPECA
All the main programs are in main_PECA.m file. Please run the script and get the result from the folder called Output. All the selection status of each RE, and TF-RE-TG triplets are listed.

# Requirements
MATLAB 2018a

# Time
It usually takes several seconds on each gene. In total, about 24 hours are required to processing all the genes, and write all text output files.
