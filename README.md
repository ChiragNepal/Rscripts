# Rscripts Repository

This directory contains R scripts and environment files for reproducible analysis in R.

---

## Environment Setup

- **`installCondaR.yml`** — Conda environment setup for R and Bioconductor packages.  
  Use this file to install the exact R version and packages required to reproduce the analysis environment.

---

## Functions

This section contains R scripts and functions applicable to multiple analyses.  
The idea is that if you provide your data objects, these scripts will generate plots and save results automatically.

- **`function_DEG_GSEA`** — Performs GSEA (GO, KEGG, HALLMARK) analysis given a DESeq2 `dds` object.
- **`plotCorrelation.R`** — Plots correlations for given data objects.

---

> More scripts and functions will be added over time.

