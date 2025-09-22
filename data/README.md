# Data for scICER
For the main clustering workflow, see [run_scICER.R](https://github.com/yerry77/scICER-workflow/blob/main/codes/run_scICER.R)

## Data source
The datasets originate from the article:  
**"scICE: enhancing clustering reliability and efficiency of scRNA-seq data with multi-cluster label consistency evaluation"**.  

All **48 datasets** described in the paper can be downloaded from Zenodo:  
[https://zenodo.org/records/15113898](https://zenodo.org/records/15113898)  

The Zenodo repository provides the data in **CSV format** (`.csv.gz`).

## Conversion to Seurat objects
To make the data more convenient for single-cell analysis in R, we provide a script that converts the downloaded `.csv.gz` files into **Seurat objects** (`.qs` format).  

The script:
- Reads the CSV files containing gene expression matrices.  
- Creates Seurat objects with expression counts and cell type annotations.  
- Saves each object as a compressed `.qs` file.  

This allows users to directly load Seurat objects for scICER analysis without reprocessing the raw CSV files.

## Example conversion script

```r
# Load required packages
library(Seurat)
library(tidyverse)
library(qs)

# Increase connection size for large CSV files
Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)

# Define input and output directories
input_dir <- "path/to/downloaded_csv_files"   # replace with your CSV download folder
output_dir <- "path/to/save_seurat_objects"   # replace with your preferred output folder
dir.create(output_dir, showWarnings = FALSE)

# List all compressed CSV files
files <- list.files(input_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)

# Process each file
for (file_path in files) {
  file_name <- basename(file_path)
  object_name <- str_remove(file_name, "_pca_dict\\.csv\\.gz$|_dict\\.csv\\.gz$|\\.csv\\.gz$")
  
  # Read raw data
  raw_data <- read_csv(file_path, col_names = FALSE)
  
  # Extract gene names and expression matrix
  gene_names <- make.unique(as.character(unlist(raw_data[1, -1])))
  expr_data <- apply(raw_data[-1, -1], 2, as.numeric)
  expr_mat <- t(expr_data)
  rownames(expr_mat) <- gene_names
  colnames(expr_mat) <- paste0("Cell", seq_len(ncol(expr_mat)))
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr_mat, project = object_name)
  
  # Add cell type annotation (if available in first column)
  celltype <- raw_data$X1[2:nrow(raw_data)]
  seurat_obj$celltype <- celltype
  
  # Save as .qs object
  qs::qsave(seurat_obj, file = file.path(output_dir, paste0(object_name, ".qs")), nthreads = 4)
}
