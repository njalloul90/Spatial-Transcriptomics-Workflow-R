#!/usr/bin/env Rscript
# scripts/04-pseudobulk-dge.R
# Pseudobulk differential expression: disease vs healthy

message("=== Running Pseudobulk DGE ===")

suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(limma)
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
# Ensure metadata has condition column
# -----------------------------
if(!"condition" %in% colnames(seurat_obj@meta.data)){
  # For example dataset, assign dummy condition
  seurat_obj$condition <- sample(c("healthy", "disease"), ncol(seurat_obj), replace = TRUE)
}

# -----------------------------
# Pseudobulk aggregation function
# -----------------------------
pseudobulk_counts <- function(seurat_obj, ct_col = "cell_type", sample_col = "orig.ident") {
  cells <- colnames(seurat_obj)
  meta <- seurat_obj@meta.data
  meta$cell <- rownames(meta)
  groups <- meta %>% mutate(group = paste0(!!sym(sample_col), "_", !!sym(ct_col)))
  group_ids <- unique(groups$group)
  mat <- seurat_obj[["RNA"]]@counts
  pb <- sapply(group_ids, function(g){
    cells_g <- groups$cell[groups$group == g]
    if(length(cells_g)==0) return(rep(0, nrow(mat)))
    rowSums(mat[, cells_g, drop=FALSE])
  })
  colnames(pb) <- group_ids
  rownames(pb) <- rownames(mat)
  list(counts = pb, meta = groups %>% distinct(group, .keep_all = TRUE))
}

pb <- pseudobulk_counts(seurat_obj)
counts <- pb$counts
col_meta <- pb$meta

# Attach condition
col_meta$condition <- seurat_obj$condition[col_meta$cell]

# Run DGE per cell type
results_list <- list()
for(ct in unique(col_meta$cell_type)){
  sel <- col_meta$group[col_meta$cell_type==ct]
  if(length(sel) < 2) next
  y <- DGEList(counts=counts[, sel, drop=FALSE])
  keep <- filterByExpr(y, group = col_meta$condition[col_meta$cell_type==ct])
  y <- y[keep,, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ 0 + col_meta$condition[col_meta$cell_type==ct])
  colnames(design) <- levels(factor(col_meta$condition[col_meta$cell_type==ct]))
  v <- voom(y, design, plot=FALSE)
  fit <- lmFit(v, design)
  if(all(c("disease","healthy") %in% colnames(design))){
    contrast <- makeContrasts(DvsH = disease - healthy, levels = design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    top <- topTable(fit2, number=Inf, sort.by="P")
    top$cell_type <- ct
    results_list[[ct]] <- top
    write.csv(top, file.path("results", paste0("DGE_", ct, ".csv")), row.names = TRUE)
  }
}

message("Pseudobulk DGE complete.")
