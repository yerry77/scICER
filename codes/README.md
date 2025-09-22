# Single-Cell Clustering Scripts

This directory contains R scripts for single-cell RNA-seq clustering analysis using scICER.

## Scripts

1. **run_scICER.R**  
   Runs scICER clustering on a Seurat object. Outputs IC plot, ECS score table, and updated Seurat object with optimal clusters.
   https://github.com/yerry77/scICER-workflow/blob/main/codes/run_scICER.R

2. **Harmony integration + scICER**  
   Performs **Harmony integration**, UMAP visualization, and scICER clustering on a Seurat object. Outputs UMAP plots, IC plot, ECS score table, and updated Seurat object with optimal clusters.
   https://github.com/yerry77/scICER-workflow/blob/main/codes/Harmony%20integration%20%2B%20scICER

4. **scVI integration + scICER**  
   Performs **scVI-based integration**, UMAP visualization, and scICER clustering on a Seurat object containing scVI latent space. Outputs UMAP plots, IC plot, ECS score table, and updated Seurat object with optimal clusters.
   https://github.com/yerry77/scICER-workflow/blob/main/codes/scVI%20integration%20%2B%20scICER

## Data information and conversion instructions
These workflows use publicly available single-cell RNA-seq datasets. Details on the data sources, preprocessing, and conversion to Seurat objects are provided in the repository’s data documentation:
https://github.com/yerry77/scICER-workflow/blob/main/data/README.md
