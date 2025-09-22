# Single-Cell Clustering Scripts

This directory contains R scripts for single-cell RNA-seq clustering analysis using scICER.

## Scripts

1. **scICER_only.R**  
   Runs scICER clustering on a Seurat object. Outputs IC plot, ECS score table, and updated Seurat object with optimal clusters.

2. **harmony_scICER.R**  
   Performs **Harmony integration**, UMAP visualization, and scICER clustering on a Seurat object. Outputs UMAP plots, IC plot, ECS score table, and updated Seurat object with optimal clusters.

3. **scvi_scICER.R**  
   Performs **scVI-based integration**, UMAP visualization, and scICER clustering on a Seurat object containing scVI latent space. Outputs UMAP plots, IC plot, ECS score table, and updated Seurat object with optimal clusters.


