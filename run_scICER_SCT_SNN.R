#!/usr/bin/env Rscript
# ===============================================================
# Script: run_scICER_SCT_SNN.R
# Purpose: Run scICER clustering on data normalized with Seurat's SCT assay (SCTransform, regressing out percent.mt) using SNN graph construction.
# Input:  qs file containing a Seurat object (pre-processed with SCT)
# Output: UMAP plots, IC plot, clustering results (tsv, qs)
# ===============================================================

library(Seurat)
library(ggplot2)
library(qs)
library(readr)
library(cowplot)
library(tidyverse)
library(mclust)
library(scICER)
library(ClustAssess)
library(uwot)

# Parse input arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide a qs file path, e.g.: Rscript run_scICER_SCT_SNN.R /path/to/seurat.scICER.qs")
}
qs_path <- args[1]
nthreads_qread <- 40

# Extract sample name from qs path (parent directory name)
sample_name <- basename(dirname(qs_path))
message("Processing sample: ", sample_name)

# Load Seurat object
seurat_obj <- qs::qread(qs_path, nthreads = nthreads_qread)
# Add mitochondrial percentage metadata
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")

# SCT normalization (regress out percent.mt)
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)

# Downstream analysis
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)

# Overwrite original qs file
qs::qsave(seurat_obj, qs_path)
message("Seurat object processed and saved (overwritten): ", qs_path)

# Output directory (example)
output_base_dir <- "results/SCT_SNN"
output_dir <- file.path(output_base_dir, sample_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save UMAP colored by cell type
pdf(file.path(output_dir, "umap.celltype.pdf"), width = 7.5, height = 6)
DimPlot(seurat_obj, group.by = "celltype")
dev.off()

# Run scICER clustering (SCT + SNN)
scice_results <- scICE_clustering(
  object = seurat_obj,
  cluster_range = 2:20,
  remove_threshold = Inf,
  n_workers = 8,
  n_trials = 30,
  n_bootstrap = 200,
  seed = 123,
  verbose = TRUE,
  graph_name = "SCT_snn"
)

# Save IC plot
ic_plot <- plot_ic(scice_results, threshold = 1.005)
ggsave(filename = file.path(output_dir, "scICER_IC_plot.pdf"),
       plot = ic_plot, device = "pdf", width = 7.5, height = 4.5)
write_tsv(ic_plot$data, file.path(output_dir, "scICER_plot_data.tsv"))

# Create dataframe with IC and ECS scores
df <- data.frame(
  cluster_number = scice_results$n_cluster,
  ic_score = scice_results$ic
)
df$is_consistent <- df$ic_score <= 1.005

sample_obj_df <- get_robust_labels(scice_results, return_seurat = FALSE, threshold = Inf)
sample_obj_df$celltype <- seurat_obj@meta.data$celltype

for (k in 2:20) {
  cluster_col <- paste0("clusters_", k)
  if (cluster_col %in% names(sample_obj_df)) {
    ari_score <- element_sim(sample_obj_df$celltype, sample_obj_df[[cluster_col]])
    df[df$cluster_number == k, "ECS_score"] <- ari_score
  } else {
    df[df$cluster_number == k, "ECS_score"] <- NA
  }
}

# Save clustering results
seurat_obj <- get_robust_labels(scice_results, return_seurat = TRUE, threshold = Inf)
write_tsv(df, file.path(output_dir, "scICER_cluster_data.tsv"))
qs::qsave(seurat_obj, file = file.path(output_dir, "seurat.scICER.qs"), nthreads = nthreads_qread)

# Generate clustering UMAP plots
dimplot_list <- list()
for (k in 2:20) {
  col_name <- paste0("clusters_", k)
  if (col_name %in% colnames(seurat_obj@meta.data)) {
    dimplot_list[[col_name]] <- DimPlot(seurat_obj, group.by = col_name, label = TRUE) +
      ggtitle(paste("Clusters:", k))
  }
}
dimplot_list[["celltype"]] <- DimPlot(seurat_obj, group.by = "celltype", label = TRUE) +
  ggtitle("celltype")

dimplot_grid <- cowplot::plot_grid(plotlist = dimplot_list, ncol = 5)
ggsave(filename = file.path(output_dir, "scICER_DimPlot_grid.pdf"),
       plot = dimplot_grid, device = "pdf", width = 32, height = 16)
