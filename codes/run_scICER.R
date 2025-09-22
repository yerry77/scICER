#!/usr/bin/env Rscript
# ===============================================================
# Input:  CSV/CSV.GZ expression matrix OR qs file OR h5ad file
# Method: RNA, SCT, or scLENS
# graph: SNN,KNN,or UMAP
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
library(SeuratDisk)
library(reticulate)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_scICER_single_sample.R <input_file> <RNA|SCT|scLENS>")
}
file_path <- args[1]
method <- toupper(args[2])  # RNA, SCT, scLENS

output_base_dir <- "/output_dir/"
nthreads_qread <- 40
sample_name <- gsub('\\.csv\\.gz$|\\.qs$|\\.h5ad$', '', basename(file_path))
message("Processing file: ", file_path, " | Method: ", method)

# Load input
if (method %in% c("RNA", "SCT")) {
  
  if (grepl("\\.qs$", file_path)) {
    seurat_obj <- qs::qread(file_path, nthreads = nthreads_qread)
  } else if (grepl("\\.csv(\\.gz)?$", file_path)) {
    expr_matrix <- read.csv(file_path, row.names = 1, check.names = FALSE)
    seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = sample_name)
  } else {
    stop("Unsupported file type for RNA/SCT. Use .csv/.csv.gz or .qs")
  }
  
  # Preprocessing
  if (method == "RNA") {
    message("Running Seurat RNA normalization pipeline...")
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)  
    seurat_obj <- ScaleData(seurat_obj)
    
  } else if (method == "SCT") {
    message("Running SCTransform pipeline (regressing out percent.mt)...")
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  }
  
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  
} else if (method == "SCLENS") {
  
  # Convert h5ad -> h5seurat
  if (!grepl("\\.h5ad$", file_path)) stop("scLENS method requires h5ad input")
  
  Convert(file_path, dest = "h5seurat", overwrite = TRUE)
  file_h5seurat <- sub("\\.h5ad$", ".h5seurat", file_path)
  seurat_obj <- LoadH5Seurat(file_h5seurat, meta.data = FALSE, misc = TRUE)
  
  # Use Python to read obs metadata
  sc <- import("scanpy")
  adata <- sc$read_h5ad(file_path)
  seurat_obj <- AddMetaData(seurat_obj, metadata = adata$obs)
  
  # PCA-based graph
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dim(seurat_obj@reductions$pca)[2], reduction = "pca")
  
} else {
  stop("Invalid method. Use RNA, SCT, or scLENS")
}

# UMAP graph construction 
pca <- Embeddings(seurat_obj, reduction = "pca")
umap_result <- uwot::umap(
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
seurat_obj[["umap"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)
rownames(umap_graph) <- rownames(seurat_obj@meta.data)
colnames(umap_graph) <- rownames(seurat_obj@meta.data)
seurat_obj[["umap_graph"]] <- as.Graph(umap_graph)


# Create output directory
output_dir <- file.path(output_base_dir, sample_name, method)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Optional plots: UMAP by celltype
if ("celltype" %in% colnames(seurat_obj@meta.data)) {
  pdf(file.path(output_dir, paste0(sample_name, "_", method, "_umap.celltype.pdf")), width = 7.5, height = 6)
  print(DimPlot(seurat_obj, group.by = "celltype"))
  dev.off()
}

# Run scICER clustering
sample_obj <- seurat_obj
scice_results <- scICE_clustering(
  object = sample_obj,
  cluster_range = 2:20,
  remove_threshold = 1.005,  # IC threshold for cluster consistency
                             # Clusters with IC > 1.005 are considered unstable and may be removed.
                             # Default 1.005 is an empirical value; you can adjust it depending on your dataset.
  n_workers = 80,
  n_trials = 15,
  n_bootstrap = 100,
  seed = 123,
  verbose = TRUE,
  graph_name = "RNA_snn"  # Graph to use for clustering:
                         # RNA method:
                         #   "RNA_snn"   - RNA method, SNN graph
                         #   "RNA_nn"    - RNA method, KNN graph
                         #   "umap_graph"  - RNA method, UMAP graph
                         # SCT method:
                         #   "SCT_snn"   - SCT method, SNN graph
                         #   "SCT_nn"    - SCT method, KNN graph
                         #   "umap_graph"- SCT method, UMAP graph
                         # scLENS method:
                         #   "RNA_snn"   - SNN graph from scLENS preprocessing
                         #   "RNA_nn"    - KNN graph from scLENS preprocessing
                         #   "umap_graph"- UMAP graph from scLENS preprocessing

)

# Save IC plot
ic_plot <- plot_ic(scice_results, threshold = 1.005)
ggsave(filename = file.path(output_dir, paste0(sample_name, "_", method, "_scICER_IC_plot.pdf")),
       plot = ic_plot, device = "pdf", width = 7.5, height = 4.5)
write_tsv(ic_plot$data, file.path(output_dir, paste0(sample_name, "_", method, "_scICER_plot_data.tsv")))

# Collect IC + ECS scores
df <- data.frame(
  cluster_number = scice_results$n_cluster,
  ic_score = scice_results$ic
)
df$is_consistent <- df$ic_score <= 1.005

seurat_obj.df <- get_robust_labels(scice_results, return_seurat = FALSE, threshold = Inf)
if ("celltype" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj.df$celltype <- seurat_obj@meta.data$celltype
  for (k in 2:20) {
    cluster_col <- paste0("clusters_", k)
    if (cluster_col %in% names(seurat_obj.df)) {
      df[df$cluster_number == k, "ECS_score"] <- element_sim(seurat_obj.df$celltype, seurat_obj.df[[cluster_col]])
    } else {
      df[df$cluster_number == k, "ECS_score"] <- NA
    }
  }
}

# Save results
sample_obj <- get_robust_labels(scice_results, return_seurat = TRUE, threshold = Inf)
write_tsv(df, file.path(output_dir, paste0(sample_name, "_", method, "_scICER_cluster_data.tsv")))
qs::qsave(sample_obj, file = file.path(output_dir, paste0(sample_name, "_", method, "_seurat.scICER.qs")), nthreads = nthreads_qread)

# Generate clustering UMAPs
dimplot_list <- list()
for (k in 2:20) {
  col_name <- paste0("clusters_", k)
  if (col_name %in% colnames(sample_obj@meta.data)) {
    dimplot_list[[col_name]] <- DimPlot(sample_obj, group.by = col_name, label = TRUE) +
      ggtitle(paste("Clusters:", k))
  }
}
if ("celltype" %in% colnames(sample_obj@meta.data)) {
  dimplot_list[["celltype"]] <- DimPlot(sample_obj, group.by = "celltype", label = TRUE) +
    ggtitle("celltype")
}
dimplot_grid <- cowplot::plot_grid(plotlist = dimplot_list, ncol = 5)
ggsave(filename = file.path(output_dir, paste0(sample_name, "_", method, "_scICER_DimPlot_grid.pdf")),
       plot = dimplot_grid, device = "pdf", width = 32, height = 16)

