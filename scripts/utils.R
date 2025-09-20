# utils.R
library(Seurat)
library(tidyverse)

# quick helper: load 10x Visium directory (Seurat)
load_visium <- function(sample_dir, sample_name){
  message("Loading: ", sample_name)
  s <- Load10X_Spatial(data.dir = sample_dir, filename = "filtered_feature_bc_matrix.h5", assay="Spatial")
  s$sample_id <- sample_name
  return(s)
}

# function to compute pseudobulk counts per sample+celltype
pseudobulk_counts <- function(seurat_obj, ct_col = "cell_type", sample_col = "sample_id"){
  require(dplyr)
  require(Matrix)
  cells <- colnames(seurat_obj)
  meta <- seurat_obj@meta.data
  stopifnot(ct_col %in% colnames(meta), sample_col %in% colnames(meta))
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
  return(list(counts = pb, meta = groups %>% distinct(group, .keep_all = TRUE)))
}
