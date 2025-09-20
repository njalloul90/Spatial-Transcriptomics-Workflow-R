#!/usr/bin/env Rscript
# scripts/05-abundance-and-regression.R
# Cell type abundance estimation & regression

message("=== Running Cell Type Abundance Analysis ===")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
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
# Compute proportions per sample
# -----------------------------
meta <- seurat_obj@meta.data %>% as.data.frame()
prop_table <- meta %>%
  group_by(orig.ident, cell_type) %>%
  tally() %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write.csv(prop_table, file.path("results", "celltype_proportions.csv"), row.names = FALSE)
message("Cell type abundance table saved.")
