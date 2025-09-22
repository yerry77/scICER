# Data for scICER

This folder provides example datasets and instructions for running **scICER**.

For the main clustering workflow, see https://github.com/yerry77/scICER-workflow/tree/main/codes.



## 1. scICE benchmark datasets (48 datasets)

### Source
These datasets are from the article:  
**"scICE: enhancing clustering reliability and efficiency of scRNA-seq data with multi-cluster label consistency evaluation"** https://www.nature.com/articles/s41467-025-60702-8#code-availability

Download from Zenodo:  
[https://zenodo.org/records/15113898](https://zenodo.org/records/15113898)

Data format: `.csv.gz`

### Conversion to Seurat objects
To facilitate R-based single-cell analysis, CSVs can be converted to **Seurat objects** (`.qs`):

```r
library(Seurat)
library(tidyverse)
library(qs)

Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)

input_dir <- "path/to/downloaded_csv_files"
output_dir <- "path/to/save_seurat_objects"
dir.create(output_dir, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)

for (file_path in files) {
  file_name <- basename(file_path)
  object_name <- str_remove(file_name, "_pca_dict\\.csv\\.gz$|_dict\\.csv\\.gz$|\\.csv\\.gz$")
  
  raw_data <- read_csv(file_path, col_names = FALSE)
  
  gene_names <- make.unique(as.character(unlist(raw_data[1, -1])))
  expr_data <- apply(raw_data[-1, -1], 2, as.numeric)
  expr_mat <- t(expr_data)
  rownames(expr_mat) <- gene_names
  colnames(expr_mat) <- paste0("Cell", seq_len(ncol(expr_mat)))
  
  seurat_obj <- CreateSeuratObject(counts = expr_mat, project = object_name)
  seurat_obj$celltype <- raw_data$X1[2:nrow(raw_data)]
  
  qs::qsave(seurat_obj, file = file.path(output_dir, paste0(object_name, ".qs")), nthreads = 40)
}
```


## 2. Integration benchmark datasets for scICER

### Source
These datasets are from the article:  
**"Benchmarking atlas-level data integration in single-cell genomics" (Nature Methods)**  
[https://www.nature.com/articles/s41592-021-01336-8#data-availability](https://www.nature.com/articles/s41592-021-01336-8#data-availability)

Download from Figshare:  
[https://figshare.com/articles/dataset/Benchmarking_atlas-level_data_integration_in_single-cell_genomics_-_integration_task_datasets_Immune_and_pancreas_/12420968](https://figshare.com/articles/dataset/Benchmarking_atlas-level_data_integration_in_single-cell_genomics_-_integration_task_datasets_Immune_and_pancreas_/12420968)

Data format: `.h5ad`

### Conversion to Seurat objects
H5AD files can be converted to **Seurat objects** (`.qs`) for convenient analysis in R:

```r
library(Seurat)
library(qs)
library(SeuratDisk)
library(reticulate)

reticulate::use_python("path/to/python")  # adjust to your Python environment
sc <- import("scanpy")

input_dir <- "path/to/downloaded_h5ad_files"
output_dir <- "path/to/save_seurat_objects"
dir.create(output_dir, showWarnings = FALSE)

files <- list.files(input_dir, pattern = "\\.h5ad$", full.names = TRUE)

for (file_path in files) {
  sample_name <- sub('\\.h5ad$', '', basename(file_path))
  message("Processing file: ", file_path)
  
  adata <- sc$read_h5ad(file_path)
  
  counts <- t(adata$layers$as_dict()[["counts"]])
  rownames(counts) <- adata$var_names$to_list()
  colnames(counts) <- adata$obs_names$to_list()
  
  seurat_obj <- CreateSeuratObject(counts = counts, meta.data = adata$obs)
  
  qs::qsave(seurat_obj, file = file.path(output_dir, paste0(sample_name, ".qs")), nthreads = 40)
}
