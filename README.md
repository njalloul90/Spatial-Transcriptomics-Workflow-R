# Spatial-Transcriptomics-Workflow-R
Complete workflow for spatial transcriptomics analysis in R


# 🧬 Spatial Transcriptomics Analysis Workflow

This repository provides a **complete, reproducible workflow in R** for analyzing **spatial transcriptomics data**.  
It includes preprocessing, clustering, cell type annotation, and differential expression analysis — with support for **cell type abundance adjustment**.

The repo is designed for:
- 🔬 Researchers working with Visium/other spatial data  
- 🛠️ Analysts who want a reproducible, containerized workflow  
- ☁️ CI/CD integration (via GitHub Actions) for automated builds and reports  

---

## ⚙️ Workflow Overview

1. **Data QC & Preprocessing**
   - Load raw spatial transcriptomics data  
   - Filter low-quality spots/cells  
   - Normalize and scale data  

2. **Clustering & Marker Detection**
   - Perform dimensionality reduction (PCA, UMAP)  
   - Cluster cells  
   - Detect marker genes per cluster  

3. **Cell Type Annotation**
   - Reference-based annotation using [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html) and `celldex` reference datasets  

4. **Differential Expression Analysis**
   - Pseudobulk aggregation per cell type  
   - Differential testing between **disease vs. healthy** using `edgeR`/`limma`  

5. **Cell Type Abundance & Confounder Regression**
   - Estimate abundance of cell types per spatial region  
   - Regress abundance covariates out of the differential expression model  

6. **Visualization & Reporting**
   - Spatial feature plots, abundance maps, DEG volcano plots  
   - Final report rendered via RMarkdown  

---

## 🚀 Quick Start

### Run with Docker (recommended)

Build the image:
```bash
docker build -t spatial-workflow .

