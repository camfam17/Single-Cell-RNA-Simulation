library(scDesign2)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(dplyr)

column <- "Donor.full" #Status, Donor.full, Admission, Ventilated
label_granularity <- "coarse" # coarse # fine # medium
gene_filter <- 0.00
zp.cutoff <- 0.8

output.folder <- "MAITSim/FitGroups6/"
scdata = readRDS("~/MAITSim/COVIDatlas_MAIT_labelled.seu.rds")

cell_grid_fold <- c(1.5, 1.5) # specify cell counts
depth_grid_fold <- c(1, 2) # specify sequencing depths

for(val in unique(scdata@meta.data[[column]])){
  cat("column", column, ": ", val, "\n")
  for(i in seq_len(length(cell_grid_fold))){
    cat("cell:", cell_grid_fold[i], " depth:", depth_grid_fold[i], "\n")

    sim.filename <- paste(output.folder, "sim_", column, "_", val, "_", cell_grid_fold[i], "cell_", depth_grid_fold[i], "depth", ".rds", sep="")
    print(sim.filename)
    if (file.exists(sim.filename)){
      next
    }
    
    rawcounts <- readRDS(paste("~/MAITSim/FitGroups6/raw_", column, val, "_label", label_granularity, "_genefilter", gene_filter, "_zpcutoff", zp.cutoff, ".rds", sep=""))
    fit.filename <- paste(output.folder, "fit_", column, val, "_label", label_granularity, "_genefilter", gene_filter, "_zpcutoff", zp.cutoff, ".rds", sep="")
    if (!file.exists(fit.filename)){
      next
    }
    fit = readRDS(fit.filename)
    
    
    n_cell_old <- sum(sapply(fit, function(x) x$n_cell))
    n_depth_old <- sum(sapply(fit, function(x) x$n_read))
    cell_count_new <- round(n_cell_old * cell_grid_fold[i])
    depth_new <- round(n_depth_old * depth_grid_fold[i])
    cat("cell old: ", n_cell_old, "  cell new :", cell_count_new)
    cat("depth old: ", n_depth_old, "  depth new:", depth_new)
    
    cell_type_prop <- prop.table(table(colnames(rawcounts)))
    cell_type_prop <- cell_type_prop[names(cell_type_prop) %in% names(fit)]
    
    print('simulating')
    sim <- simulate_count_scDesign2(
      model_params   = fit,
      n_cell_new     = cell_count_new,
      cell_type_prop = cell_type_prop,
      total_count_new = depth_new,
      sim_method     = "copula",
      reseq_method   = "mean_scale", 
      cell_sample    = FALSE
    )
    rownames(sim) <- rownames(rawcounts)
    simcount1 <- Matrix(sim, sparse = TRUE)


    saveRDS(sim, file = sim.filename)
    print("data saved")
    
    
  }
  
}