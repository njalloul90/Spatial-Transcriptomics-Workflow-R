#!/usr/bin/env Rscript
# scripts/00-setup.R
# Setup script for Spatial Transcriptomics workflow
# - Installs required packages
# - Downloads small example data (stxBrain) if no user data provided
# - Ensures folder structure exists

message("=== Setting up environment for Spatial Transcriptomics workflow ===")

# -----------------------------
# Package installation function
# -----------------------------
install_if_missing <- function(pkgs, bioc = FALSE) {
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkgs, ask = FALSE, update = TRUE)
  } else {
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing) > 0) {
      install.packages(missing, repos = "https://cloud.r-project.org")
    }
  }
}

# -----------------------------
# Install required packages
# -----------------------------
cran_pkgs <- c(
  "here", "Matrix", "cowplot", "patchwork",
  "jsonlite", "scales", "rmarkdown", "SeuratData"
)

bioc_pkgs <- c(
  "Seurat", "SeuratDisk",
  "SingleCellExperiment", "SingleR", "celldex",
  "edgeR", "limma", "SpatialExperiment"
)

message("Installing CRAN packages...")
install_if_missing(cran_pkgs, bioc = FALSE)

message("Installing Bioconductor packages...")
install_if_missing(bioc_pkgs, bioc = TRUE)

# Load packages
sapply(c(cran_pkgs, bioc_pkgs), require, character.only = TRUE)

# -----------------------------
# Ensure folder structure
# -----------------------------
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("results")) dir.create("results")

# -----------------------------
# Example dataset setup
# -----------------------------
example_file <- "data/example_spatial.rds"

if (!file.exists(example_file)) {
  message("No user data found in data/. Downloading example dataset (stxBrain)...")
  SeuratData::InstallData("stxBrain")
  data("stxBrain", package = "stxBrain")
  saveRDS(stxBrain, file = example_file)
  message("Example dataset saved to ", example_file)
} else {
  message("Example dataset already exists at ", example_file)
}

# -----------------------------
# Data availability check
# -----------------------------
user_files <- list.files("data", pattern = "\\.rds$", full.names = TRUE)
if (length(user_files) > 1) {
  message("User data detected in /data: ", paste(basename(user_files), collapse = ", "))
  message("Workflow will use user data instead of example dataset.")
} else {
  message("Using example dataset. To analyze your own data, place .rds files in /data.")
}

message("=== Setup complete! ===")
