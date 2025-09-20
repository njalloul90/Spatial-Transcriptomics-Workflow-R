#!/usr/bin/env Rscript
# scripts/03-celltype-annotation.R
# Cell type annotation using SingleR

message("=== Running Cell Type Annotation ===")

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
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
# SingleR annotation
# -----------------------------
message("Converting to SingleCellExperiment...")
sce <- as.SingleCellExperiment(seurat_obj, assay = "SCT")
ref <- celldex::HumanPrimaryCellAtlasData()
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
seurat_obj$SingleR_label <- pred$labels

# Optional: propagate cluster-level labels
seurat_obj$cell_type <- seurat_obj$SingleR_label

saveRDS(seurat_obj, file.path("results", "seurat_annotated.rds"))
write.csv(table(seurat_obj$cell_type, seurat_obj$seurat_clusters),
          file.path("results", "celltype_counts.csv"))
message("Cell type annotation complete.")
