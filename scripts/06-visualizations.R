#!/usr/bin/env Rscript
# scripts/06-visualizations.R
# UMAP and spatial visualizations

message("=== Generating Visualizations ===")

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(cowplot)
  library(ggplot2)
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
# UMAP
# -----------------------------
if("umap" %in% names(seurat_obj@reductions)){
  p1 <- DimPlot(seurat_obj, reduction="umap", group.by="cell_type", label=TRUE) + ggtitle("UMAP Cell Types")
  ggsave(file.path("results", "umap_celltypes.png"), p1, width=8, height=6)
}

# -----------------------------
# Spatial plot (example feature)
# -----------------------------
feature <- rownames(seurat_obj)[1]  # top gene
p2 <- SpatialFeaturePlot(seurat_obj, features=feature) + ggtitle(feature)
ggsave(file.path("results", paste0("spatial_", feature, ".png")), p2, width=8, height=6)

message("Visualizations saved to results/")
