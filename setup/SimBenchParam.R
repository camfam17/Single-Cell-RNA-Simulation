
library(scDesign2)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(plyr)
library(dplyr)
library(Seurat)
library(sctransform)
library(igraph)
library(factoextra)
library(ComplexHeatmap)
library(circlize)
library(EpicTools)
require(Hmisc)
require(dplyr)
require(openxlsx)
require(ggplot2)
library(ggpubr)
require(cowplot)
library(data.table)
library(RColorBrewer)
library(rowr)
library(SingleR)
library(scater)
library(pheatmap)
library(nichenetr)
library(tidyverse)

library(Matrix)
library(Matrix.utils)
library(Seurat)
library(dplyr)
library(SimBench)
library(parallel)
library(DESeq2)



target <- "All" # COVID Healthy
column <- "Status"
# target <- "H1" # C2
# column <- "Donor"
label_granularity <- "fine" # coarse # fine
gene_filter <- 0.00
zp.cutoff <- 0.8



print("Loading data")

# realmat <- readRDS(file=paste("~/MAITSim/rawmats/raw_", column, target, "_label", label_granularity, "_genefilter", gene_filter, "_zpcutoff", zp.cutoff, ".rds", sep=""))
realmat <- readRDS(file=paste("~/MAITSim/out4/" , "1to1raw_", column, "_", target, ".rds", sep=""))


real_labels <- as.character(colnames(realmat))
colnames(realmat) <- make.unique(real_labels, sep = "_cell")

# realmat <- as.matrix(realmat)
realsce <- SingleCellExperiment(assays = list(counts = realmat))

print(realmat[1:10, 1:10])


# simmat <- readRDS(file=paste("~/MAITSim/simmats/data_", column, target, "_lbl", label_granularity, "_gf", gene_filter, "_zp", zp.cutoff, ".rds", sep=""))
simmat <- readRDS(file=paste("~/MAITSim/out4/", "1to1sim_", column, "_", target, ".rds", sep=""))

sim_labels  <- as.character(colnames(simmat))
colnames(simmat)  <- make.unique(sim_labels,  sep = "_cell")

# simmat <- as.matrix(simmat)
simsce <- SingleCellExperiment(assays = list(counts = simmat))

print(simmat[1:10, 1:10])





print("Evaluating Parameters")


colData(realsce)$celltype <- factor(real_labels)
colData(simsce)$celltype  <- factor(sim_labels)

parameter_result <- eval_parameter(real = realsce, sim = simsce, type = "raw" , method = "samplemethod")


print("Saving RDS")

# saveRDS(parameter_result, paste("~/MAITSim/SimBench/SimBenchParameter", column, target, "_lbl", label_granularity, ".rds"))
saveRDS(parameter_result, paste("~/MAITSim/out4/SimBenchParameter_", column, "_", target, ".rds"))


print("Saved RDS")
print("FINISHED")
