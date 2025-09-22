# scICER: A Flexible R Package for Assessing Clustering Stability in Integrated Single-Cell Analyses

## Introduction
**scICER** is an R package for evaluating cluster consistency in single-cell RNA-seq data. It reimplements the scICE methodology to assess the stability of cluster labels across stochastic runs. Fully integrated with Seurat, scICER works on embeddings corrected by Harmony or scVI, enabling consistent clustering analysis in complex, multi-sample datasets.

## data sources, preprocessing, and conversion
This folder provides example datasets and instructions for running scICER.
https://github.com/yerry77/scICER-workflow/blob/main/data/README.md

## codes examples to run scICER 
### This folder contains example R scripts demonstrating the use of **scICER**, a tool for consistent clustering of single-cell RNA-seq data. 
The scripts include a basic usage example as well as two workflows showing how scICER can be combined with common integration methods (**Harmony** and **scVI**) for downstream analysis. Each script produces IC plots, ECS scores, and updates the Seurat object with optimal clusters. https://github.com/yerry77/scICER/tree/main/codes
