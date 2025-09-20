# 00-setup.R
# Install / load packages, set file paths

required_pkgs <- c(
  "Seurat", "SeuratDisk", "SingleCellExperiment", "SingleR", "celldex",
  "edgeR", "limma", "scales", "Matrix", "tidyverse", "patchwork",
  "cowplot", "jsonlite"
)

install_if_missing <- function(pkgs){
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if(length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
}
install_if_missing(required_pkgs)

# Bioconductor packages
if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("SingleR", "celldex", "SpatialExperiment")
for(p in bioc_pkgs) if(!p %in% installed.packages()[,1]) BiocManager::install(p, ask = FALSE)

# load
library(Seurat)
library(SeuratDisk)
library(SingleR)
library(celldex)
library(edgeR)
library(limma)
library(tidyverse)
library(patchwork)
library(cowplot)

# Paths - edit to your environment
root_dir <- here::here()
data_dir <- file.path(root_dir, "data")
raw_dir <- file.path(data_dir, "raw")
proc_dir <- file.path(data_dir, "processed")
results_dir <- file.path(root_dir, "results")
fig_dir <- file.path(results_dir, "figures")
tables_dir <- file.path(results_dir, "tables")
dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

message("Setup complete. Edit paths in 00-setup.R if needed.")
