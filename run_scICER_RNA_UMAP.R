#!/usr/bin/env Rscript
# ===============================================================
# Script: run_scICER_RNA_UMAP.R
# Purpose: Run scICER clustering on data normalized with Seurat's RNA assay (NormalizeData) using UMAP graph for clustering.
# Input:  CSV/CSV.GZ expression matrix (genes × cells) OR qs file containing Seurat object.
# Output: UMAP plots, IC plot, clustering results (tsv, qs).
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
library(igraph)

Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)

# Parse input arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_scICER_RNA_UMAP.R <input.csv.gz|input.qs>")
}
file_path <- args[1]

# Output settings
output_base_dir <- "results_RNA_UMAP"
nthreads_qread <- 40
sample_name <- gsub('\\.csv\\.gz$|_pca_dict\\.csv\\.gz|\\.qs$', '', basename(file_path))
message("Processing file: ", file_path)

# Load Seurat object
seurat_filtered <- qs::qread(file_path, nthreads = nthreads_qread)

# Standard RNA normalization pipeline
seurat_filtered <- NormalizeData(seurat_filtered)
seurat_filtered <- FindVariableFeatures(seurat_filtered)  
seurat_filtered <- ScaleData(seurat_filtered)
seurat_filtered <- RunPCA(seurat_filtered)

# Create UMAP coordinates and UMAP-based graph
seurat_obj <- seurat_filtered
pca <- Embeddings(seurat_obj, reduction = "pca")

umap_result <- umap(
  pca,
  ret_model = TRUE,
  ret_extra = c("fgraph"),
  min_dist = 0.1,
  metric = "cosine",
  n_neighbors = 15L,
  n_components = 2L
)

umap_coords <- umap_result$embedding
umap_graph <- umap_result$fgraph

# Add UMAP coordinates to Seurat object
seurat_obj[["umap"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)

# Convert fgraph to igraph and assign as "umap_graph"
row.names(umap_graph) <- row.names(seurat_obj@meta.data)
colnames(umap_graph) <- row.names(seurat_obj@meta.data)
seurat_obj[["umap_graph"]] <- as.Graph(umap_graph)

# Output directory
output_dir <- file.path(output_base_dir, sample_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save UMAP colored by cell type
pdf(file.path(output_dir, "umap.celltype.pdf"), width = 7.5, height = 6)
DimPlot(seurat_filtered, group.by = "celltype")
dev.off()

# Run scICER clustering
sample_obj <- seurat_obj
scice_results <- scICE_clustering(
  object = sample_obj,
  cluster_range = 2:20,
  remove_threshold = Inf,
  n_workers = 8,
  n_trials = 15,
  n_bootstrap = 100,
  seed = 123,
  verbose = TRUE,
  graph_name = "umap_graph"
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

# Evaluate ECS scores
seurat_filtered.df <- get_robust_labels(scice_results, return_seurat = FALSE, threshold = Inf)
seurat_filtered.df$celltype <- seurat_filtered@meta.data$celltype

for (k in 2:20) {
  cluster_col <- paste0("clusters_", k)
  if (cluster_col %in% names(seurat_filtered.df)) {
    ari_score <- element_sim(seurat_filtered.df$celltype, seurat_filtered.df[[cluster_col]])
    df[df$cluster_number == k, "ECS_score"] <- ari_score
  } else {
    df[df$cluster_number == k, "ECS_score"] <- NA
  }
}

# Save results
sample_obj <- get_robust_labels(scice_results, return_seurat = TRUE, threshold = Inf)
write_tsv(df, file.path(output_dir, "scICER_cluster_data.tsv"))
qs::qsave(sample_obj, file = file.path(output_dir, "seurat.scICER.qs"), nthreads = nthreads_qread)

# Generate clustering UMAPs
dimplot_list <- list()
for (k in 2:20) {
  col_name <- paste0("clusters_", k)
  if (col_name %in% colnames(sample_obj@meta.data)) {
    dimplot_list[[col_name]] <- DimPlot(sample_obj, group.by = col_name, label = TRUE) +
      ggtitle(paste("Clusters:", k))
  }
}
dimplot_list[["celltype"]] <- DimPlot(sample_obj, group.by = "celltype", label = TRUE) +
  ggtitle("celltype")

dimplot_grid <- cowplot::plot_grid(plotlist = dimplot_list, ncol = 5)
ggsave(filename = file.path(output_dir, "scICER_DimPlot_grid.pdf"),
       plot = dimplot_grid, device = "pdf", width = 32, height = 16)
