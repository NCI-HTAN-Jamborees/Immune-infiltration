# Immune-infiltration

### Human Tumor Atlas Network (HTAN) Data Jamboree | Dec. 4-7, 2023 

## Overview

Immune infiltration is a higly complex and orchestrated process in tumor initation and progression. The aim of the study is to utilize the paired transcriptomic and spatial transcriptomic data available through the HTAN initiative to quantify the relative infiltration of immune cells within the tumors and discover relevant biological insights. 

We focused on utilizing data from patients with matched spatial transcriptomic (10X Visium platform) and single-cell transcriptomic (scRNA-seq) profiling in order to assay not only the immune infiltration (primarily through spatial neighborhood analysis) but also to uncover biological programs and crosstalk underlying immune infiltration as well. 

We utilized a larger scRNA-seq dataset (n=137,000 cells pre QC) to label coarse (broad) and fine (specific) immune cell types. We applied Cottrazm pipeline for spatial transcriptomics datasets to determine putative malignant, boundary and non malignant regions. This is accomplished using InferCNV tool, and integrated with Seurat clustering. Immune cells were used as a normal reference to determine malignant spots. After this Seuerat spot-based scRNA-seq reference-guided deconvolution was performed.

## Workflow 

1. Data compiling and pre-processing
2. Cottrazm pipeline for morphological expression adjustment and boundary detection, single sample
3. Seurat pipeline for scRNA-seq annotation of large cohort
4. Seurat pipeline for combined scRNA-seq/Visium analysis
5. Deconvolution of Visium spots using scRNA-seq annotated data
6. Evaluation of immune infiltration

## Future directions

## Description of code utilized:

`HTAN_Wrangling.Rmd`: R Markdown code to analyze the metadata from the HTAN portal and identify which samples / studies have paired spatial/single-cell transcriptomic data from patients. We identify 13 patients from HTAN Vanderbilt dataset and 10 patients from HTAN WUSTL dataset and match them with relevant file names for downstream analyses.

`metadata_analysis.Rmd`: R Markdown code to analyze the metadata from the HTAN portal and use the immune subset enrichments identified from the other analyses to measure patient-level outcomes, such as overall heterogeneity and survival (specifically, progression-free survival) outcomes.

spatial/cottrazm/ folder: 
* `conda-instr.txt` - instructions on how to install necessary conda environment into SevenBridges DataStudio
* `cottrazm_pipeline.Rmd` - R Markdown code to run the Cottrazm pipeline
* `install.R` - installation requirements, run before pipeline
* `read_custom.py` - modified SME normalization script that is run during morphological normalization step in cottrazm_pipeline.Rmd

scRNA-seq folder:
* `gene_module_score.Rmd` - R markdown code to get the gene module score plots, essential input and output files are on CGC Seven Bridge portal

## Team

- **Noah Lee**,  Genentech, South San Francisco CA - gCS and Translational Medicine | Oncology
- **Jason Weirather**, Dana-Farber Cancer Institute, Boston MA - Department of Data Science Center for Immuno-Oncology
- **Rena Elkin**, Memorial Sloan Kettering Cancer Center - Medical Physics Department
- **Hope Mumme**, Emory University, Atlanta GA - Emory School of Medicine Department of Biomedical Informatics
- **Archana Balan**, Johns Hopkins University - Biomedical Engineering Program
- **Nina Steele**, Henry Ford Hospital, Detroit MI - Department of Surgery, Henry Ford Pancreatic Cancer Center
- **Lucy Zhengyi Chen**, Vanderbilt University, Nashville TN - Chemical and Physical Biology Program
