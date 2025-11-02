library(scDesign2)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(dplyr)

label_map1 <- c(
  "RBC" = "RBC",
  "B" = "B",
  "PB" = "PB",
  "CD14 Monocyte" = "Monocyte",
  "CD16 Monocyte" = "Monocyte",
  "CD8 T" = "T_all",
  "CD4 T" = "T_all",
  "gd T" = "T_all",
  "Platelet" = "Platelet",
  "NK" = "NK",
  "Granulocyte" = "Granulocyte",
  "DC" = "DC",
  "pDC" = "DC"
)

coarsen_labels <- function(x) {
  out <- unname(label_map1[x])     # vectorized lookup
  out[is.na(out)] <- x[is.na(out)] # keep originals if not in map
  out
}

print("Loading data")
covid_combined.nc = readRDS("~/RFolder2/downloaded.blish_covid.seu.rds")

print("Organising cell labels")
# Rename and add coarsened labels
covid_combined.nc$cell.type <- NULL
colnames(covid_combined.nc@meta.data)[colnames(covid_combined.nc@meta.data) == "cell.type.coarse"] <- "cell.type.medium"
covid_combined.nc$cell.type.coarse <- coarsen_labels(covid_combined.nc$cell.type.medium)



target <- "COVID" # COVID
column <- "Status"
# target <- "H1" # C2
# column <- "Donor"
label_granularity <- "medium" # coarse # fine # medium
gene_filter <- 0.00
zp.cutoff <- 0.8



print("Subsetting Suerat obj")
if(tolower(target) != 'all'){
  keep_cells <- covid_combined.nc@meta.data[[column]] == target
  covid_combined.nc <- subset(covid_combined.nc, cells = colnames(covid_combined.nc)[keep_cells])
}


print("Extracting and labeling count matrix")
# Extract count matrix from Seurat object
rawcounts <- GetAssayData(covid_combined.nc, assay = "RNA", slot = "counts")

# Label count matrix with cell type for scDesign2
labels <- covid_combined.nc@meta.data[colnames(rawcounts), paste("cell.type", label_granularity, sep=".")]
colnames(rawcounts) <- labels


print("Filtering low genes")
# Filter low count genes
keep <- Matrix::rowSums(rawcounts > 0) >= gene_filter * ncol(rawcounts)  # Remove genes that are expressed in less than 1% of cells
rawcounts <- rawcounts[keep, , drop=FALSE]

print("Saving extracted rawcount")
saveRDS(rawcounts, file = paste("~/RFolder2/scDesign2ExtractedRawcounts/raw_", column, target, "_label", label_granularity, "_genefilter", gene_filter, "_zpcutoff", zp.cutoff, ".rds", sep=""))


set.seed(1)
RNGkind("L'Ecuyer-CMRG")   # good practice for parallel RNG
print("fitting model")
fit <- fit_model_scDesign2(data_mat = rawcounts,
                           cell_type_sel = sort(unique(colnames(rawcounts))),
                           sim_method = "copula",
                           marginal = "auto_choose",
                           zp_cutoff = zp.cutoff,
                           ncores = min(5, length(unique(colnames(rawcounts))))
)
print("model fitted")
saveRDS(fit, file = paste("~/RFolder2/scDesign2FitModels/fit_", column, target, "_label", label_granularity, "_genefilter", gene_filter, "_zpcutoff", zp.cutoff, ".rds", sep=""))
print("model saved")
