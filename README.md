# Immune-infiltration

Immune infiltration is a higly complex and orchestrated process in tumor initation and progression. The aim of the study is to utilize the paired transcriptomic and spatial transcriptomic data available through the HTAN initiative to quantify the relative infiltration of immune cells within the tumors and discover relevant biological insights. 

We focused on utilizing data from patients with matched spatial transcriptomic (10X Visium platform) and single-cell transcriptomic (scRNA-seq) profiling in order to assay not only the immune infiltration (primarily through spatial neighborhood analysis) but also to uncover biological programs and crosstalk underlying immune infiltration as well. 

We utilized a larger scRNA-seq dataset (n=137,000 cells pre QC) to label coarse (broad) and fine (specific) immune cell types. We applied Cottrazm pipeline for spatial transcriptomics datasets to determine putative malignant, boundary and non malignant regions. This is accomplished using InferCNV tool, and integrated with Seurat clustering. Immune cells were used as a normal reference to determine malignant spots. After this Seuerat spot-based scRNA-seq reference-guided deconvolution was performed.

Description of code utilized:

HTAN_Wrangling.Rmd: R Markdown code to analyze the metadata from the HTAN portal and identify which samples / studies have paired spatial/single-cell transcriptomic data from patients. We identify 13 patients from HTAN Vanderbilt dataset and 10 patients from HTAN WUSTL dataset and match them with relevant file names for downstream analyses.

metadata_analysis.Rmd: R Markdown code to analyze the metadata from the HTAN portal and use the immune subset enrichments identified from the other analyses to measure patient-level outcomes, such as overall heterogeneity and survival (specifically, progression-free survival) outcomes.

