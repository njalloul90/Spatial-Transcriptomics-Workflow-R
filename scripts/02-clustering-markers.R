#!/usr/bin/env Rscript
# scripts/02-clustering-markers.R
# Clustering and marker detection

message("=== Running Clustering & Marker Detection ===")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(here)
})

# -----------------------------
# Load data
# -----------------------------
data_dir <- "data"
example_file <- file.path(data_dir, "example_spatial.rds")
user_files <- setdiff(list.files(data_dir, pattern = "\\.rds$", full.names = TRUE), example_file)
input_file <- if(length(user_files) > 0) user_files[1] else example_file
message("Loading dataset: ", basename(input_file))

seurat_obj <- readRDS(input_file)

# -----------------------------
# PCA, clustering, UMAP
# -----------------------------
message("Running PCA, neighbors, clustering, UMAP...")
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# -----------------------------
# Marker detection
# -----------------------------
message("Finding cluster markers...")
Idents(seurat_obj) <- seurat_obj$seurat_clusters
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers, file.path("results", "cluster_markers.csv"), row.names = FALSE)

saveRDS(seurat_obj, file.path("results", "seurat_clusters.rds"))
message("Clustering and markers saved.")
