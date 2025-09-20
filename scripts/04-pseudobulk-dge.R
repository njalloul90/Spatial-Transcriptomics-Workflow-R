# 04-pseudobulk-dge.R
# Differential expression: disease vs healthy per cell type using pseudobulk + edgeR/limma-voom

source("scripts/00-setup.R")
source("scripts/utils.R")
seurat_annot <- readRDS(file.path(proc_dir, "seurat_merged_annotated.rds"))

# metadata should contain sample_id and condition (disease/healthy)
if(!"condition" %in% colnames(seurat_annot@meta.data)){
  stop("seurat object must have a 'condition' column in metadata (e.g. 'disease' or 'healthy')")
}

# pseudobulk - create counts by sample x cell_type
pb <- pseudobulk_counts(seurat_annot, ct_col = "cell_type", sample_col = "sample_id")
counts <- pb$counts # genes x pseudobulk-samples
col_meta <- pb$meta
# parse sample and celltype from group
col_meta <- col_meta %>% mutate(
  sample = stringr::str_extract(group, "^[^_]+"),
  cell_type = stringr::str_extract(group, "[^_]+$") # after last underscore; adjust if underscores in names
)
# attach condition
# make a sample->condition mapping from original meta
sample_condition <- seurat_annot@meta.data %>% distinct(sample_id, condition) %>% rename(sample = sample_id)
col_meta <- left_join(col_meta, sample_condition, by = "sample")
# ensure columns order match counts
stopifnot(all(col_meta$group == colnames(counts)))

# For each cell type, run edgeR/limma
results_list <- list()
celltypes <- unique(col_meta$cell_type)
for(ct in celltypes){
  message("Processing cell type: ", ct)
  sel <- col_meta$group[col_meta$cell_type == ct]
  if(length(sel) < 2){
    message("  skipping ", ct, " (not enough pseudobulk samples)")
    next
  }
  y <- DGEList(counts = counts[, sel, drop=FALSE])
  keep <- filterByExpr(y, group = col_meta$condition[col_meta$cell_type == ct])
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ 0 + condition, data = col_meta[col_meta$cell_type==ct, ])
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  # contrast disease - healthy
  if(!all(c("conditiondisease","conditionhealthy") %in% colnames(design))){
    # safer: draw levels
    cond <- factor(col_meta$condition[col_meta$cell_type==ct])
    design <- model.matrix(~0 + cond)
    colnames(design) <- levels(cond)
  }
  contrast <- makeContrasts(DvsH = disease - healthy, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  top <- topTable(fit2, number = Inf, sort.by = "P")
  top$cell_type <- ct
  results_list[[ct]] <- top
  write.csv(top, file.path(tables_dir, paste0("DGE_", ct, ".csv")), row.names = TRUE)
}

# combine results table
all_res <- bind_rows(lapply(results_list, function(x) { x %>% rownames_to_column("gene") }))
write.csv(all_res, file.path(tables_dir, "DGE_all_celltypes.csv"), row.names = FALSE)
