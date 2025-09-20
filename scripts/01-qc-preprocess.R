#!/usr/bin/env Rscript
# scripts/01-qc-preprocess.R
# QC and preprocessing for Spatial Transcriptomics workflow

message("=== Running QC & Preprocessing ===")

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(Matrix)
})

# -----------------------------
# Load data (user or example)
# -----------------------------
data_dir <- "data"
example_file <- file.path(data_dir, "example_spatial.rds")

user_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
user_files <- setdiff(user_files, example_file)  # exclude example if user data present

if (length(user_files) > 0) {
  input_file <- user_files[1]
  message("User dataset found: ", basename(input_file))
} else {
  input_file <- example_file
  message("No user dataset found. Using example dataset: ", basename(input_file))
}

seurat_obj <- readRDS(input_file)

# -----------------------------
# QC filtering
# -----------------------------
message("Performing QC filtering...")

# Example QC thresholds (adjust as needed for your own data)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

qc_metrics <- data.frame(
  nCount_RNA = seurat_obj$nCount_RNA,
  nFeature_RNA = seurat_obj$nFeature_RNA,
  percent.mt = seurat_obj$percent.mt
)

# Filter cells (basic thresholds)
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           percent.mt < 15
)

# -----------------------------
# Normalization & scaling
# -----------------------------
message("Normalizing and scaling data...")
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)

# Save preprocessed object
out_file <- file.path("results", "preprocessed_seurat.rds")
saveRDS(seurat_obj, out_file)
message("Preprocessed data saved to ", out_file)

message("=== QC & Preprocessing complete! ===")
