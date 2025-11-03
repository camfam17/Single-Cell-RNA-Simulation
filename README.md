# Simulated single-cell RNA sequencing datasets to inform study design for rare cell type experiments  
  
### This project uses:  
 - the COVID-19 single-cell atlas from https://pmc.ncbi.nlm.nih.gov/articles/PMC7382903/   
 - and the scDesign2 simulator from https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02367-2  
 - to simulate single-cell RNA at different cell counts and sequencing depths to analyse how those parameters affect the detectability and downstream analysis of Mucosal-Associated Invariant T (MAIT) cells  
  
  
  
## Analysis Notebooks:
**Part 1** - Real vs Simulated Data: Notebook evaluating the fidelity of simulated scRNA-seq data against real data using SimBench, including structural, quantitative, biological, and gene-signature fidelity analyses.  
**Part 2** - Grid Search: Notebook performing a grid search over cell count and sequencing depth to assess simulation design effects on data fidelity and rare cell (MAIT) detectability.  
**Part 3** - Application/Replication: Notebook replicating downstream MAIT cell analyses (gene signatures, differential expression, gene panels) on real and simulated datasets to assess biological interpretability.    
  
Each notebook has been knit into an HTML file, so you can see the executed code and its outputs without running the code yourself.  
Executing the code with the full COVID-19 dataset can require up to 60GB of memory. The data can be subset to, for example, a single patient (Donor == C1) to allow for much faster and leaner execution - though the results will obviously differ from the main study.  
e.g.:  
```{r}
target <- "H1"
column <- "Donor"
keep_cells <- covid_combined.nc@meta.data[[column]] == target
covid_combined.nc <- subset(covid_combined.nc, cells = colnames(covid_combined.nc)[keep_cells])
```

  
## Data:  
The first step is to run the setup/LabelMAIT.Rmd notebook to download the COVID-19 single cell atlas Seurat object and organise the cell type annotations for the remainder of the analyses.  
The rest of the data, including trained simulator models, simulated count matrices and other time consuming analyses are accesable here: https://drive.google.com/file/d/1szpmLp_7kbIf9w-rZAg75QrC24xcWZ27  
This data folder contains intermediate data files - outputs of time-consuming computational runs (such as trained models, dataframes) that have been saved to disk and can be loaded by the notebooks instead of running those functions yourself. This was done to allow the code to be executed in a relatively short timeframe without having to train models yourself.  
e.g.:  
```{r}
if (file.exists(fit.filename)){
  fit <- readRDS(fit.filename)
}else{
  fit the model to the data...  
}
```


You are welcome to hide those files so that they will be generated from scratch. Please note in some cases, like training the scDesign2 model on the full dataset, can take 12-24 hours.  Simulations can take up to 20 minutes on the full dataset at the higher sequencing depths.  
The files in the "setup/" folder contain R scripts that can be executed in a non-interactive setting (unlike notebooks). These scripts can be use to train and fit models, run certain analyses and save them to disk to later be loaded by the notebooks.  

  


## Figures  
Figures that appear in the report are in "data/plots"   
