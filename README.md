Simulated single-cell RNA sequencing datasets to inform study design for rare cell type experiments

This project uses the COVID-19 single-cell atlas from https://pmc.ncbi.nlm.nih.gov/articles/PMC7382903/ 
and the scDesign2 simulator from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02367-2
to simulate single-cell RNA at different cell counts and sequencing depths to analyse how those parameters affect the detectability and downstream analysis of Mucosal-Associated Invariant T (MAIT) cells


Analysis Notebooks:


# Data:
Before doing anything, you need to run the scripts/LabelMAIT.R script to download the Seurat object and organise the cell type annotations for the remainder of this project.
The rest of the data, including trained/fit simulator models, simulated count matrices and other time consuming analyses are accesable here: 

