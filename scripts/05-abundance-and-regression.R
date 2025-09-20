# 05-abundance-and-regression.R
# Estimate cell-type abundance per region/sample, include as covariate or regress out

source("scripts/00-setup.R")
source("scripts/utils.R")
seurat_annot <- readRDS(file.path(proc_dir, "seurat_merged_annotated.rds"))

# 1) compute cell-type proportions per sample (or per region if your metadata has region)
meta <- seurat_annot@meta.data %>% as.data.frame()
prop_table <- meta %>%
  group_by(sample_id, cell_type) %>%
  tally() %>%
  group_by(sample_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write.csv(prop_table, file.path(tables_dir, "celltype_proportions_by_sample.csv"), row.names = FALSE)

# Example approach A: include cell-type proportions as covariates in the pseudobulk model
# Build a wide matrix of proportions (sample x celltype)
prop_wide <- prop_table %>% select(sample_id, cell_type, prop) %>%
  pivot_wider(names_from = cell_type, values_from = prop, values_fill = 0)
write.csv(prop_wide, file.path(tables_dir, "celltype_proportions_wide.csv"), row.names = FALSE)

# For pseudobulk DGE per cell type, include proportion of other cell types in the design:
# Re-run pseudobulk for one cell type example with covariates
pb <- pseudobulk_counts(seurat_annot, ct_col = "cell_type", sample_col = "sample_id")
counts <- pb$counts
col_meta <- pb$meta %>% mutate(
  sample = stringr::str_extract(group, "^[^_]+"),
  cell_type = stringr::str_extract(group, "[^_]+$")
)
sample_condition <- seurat_annot@meta.data %>% distinct(sample_id, condition) %>% rename(sample = sample_id)
col_meta <- left_join(col_meta, sample_condition, by = "sample")

# choose celltype to illustrate
ct_example <- unique(col_meta$cell_type)[1]
sel_groups <- col_meta$group[col_meta$cell_type == ct_example]
counts_sel <- counts[, sel_groups, drop=FALSE]
meta_sel <- col_meta[col_meta$cell_type == ct_example, ]

# attach proportion covariates (other cell types) based on sample-level prop_wide
prop_wide <- prop_wide %>% rename(sample = sample_id)
meta_sel <- left_join(meta_sel, prop_wide, by = "sample")

# prepare design: condition + other cell-type proportions (excluding the focal cell type to avoid collinearity)
focal_col <- ct_example
prop_covariates <- setdiff(colnames(prop_wide), "sample")
prop_covariates <- setdiff(prop_covariates, focal_col)
# build formula string
covars <- paste(c("condition", prop_covariates), collapse = " + ")
design <- model.matrix(as.formula(paste("~ 0 +", covars)), data = meta_sel)

# run voom/limma as before
y <- DGEList(counts = counts_sel)
keep <- filterByExpr(y, group = meta_sel$condition)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
v <- voom(y, design)
fit <- lmFit(v, design)
# contrast for condition effect (disease - healthy)
# create contrast vector dynamically; ensure names match
if(!("conditiondisease" %in% colnames(design))) {
  # fallback: find columns with prefix condition
  cond_cols <- grep("^condition", colnames(design), value = TRUE)
  if(length(cond_cols) >= 2){
    contrast <- makeContrasts(contrasts = paste0(cond_cols[1], "-", cond_cols[2]), levels = design)
  } else stop("Condition columns not found in design")
} else {
  contrast <- makeContrasts(conditiondisease - conditionhealthy, levels = design)
}
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
res <- topTable(fit2, number = Inf)
write.csv(res, file.path(tables_dir, paste0("DGE_with_prop_covariates_", focal_col, ".csv")))
message("Saved DGE with abundance covariates for ", focal_col)

# Example approach B: residualize (regress out) abundance effects per gene before DGE
# construct linear model of expression across pseudobulk samples to remove effect of proportions
# caution: residualizing expression should be done carefully; preserve biological signal
# Here demonstrate linear model residualization for clarity
expr_mat <- as.matrix(v$E) # voom-normalized logCPM
# regress out proportion covariates for each gene
X_prop <- as.matrix(meta_sel[, prop_covariates, drop=FALSE])
X_prop <- cbind(Intercept = 1, X_prop)
resid_expr <- apply(expr_mat, 1, function(gene){
  coefs <- lm.fit(X_prop, gene)$coefficients
  fitted <- X_prop %*% coefs
  gene - fitted
})
resid_expr <- t(resid_expr)
# now use resid_expr for DGE comparing condition (simple design: condition)
design_simple <- model.matrix(~0 + meta_sel$condition)
v2 <- voomList(list(Dummy=DGEList(counts=counts_sel)))$Dummy # dummy to reuse structure
# NOTE: we'll create a fake E matrix for limma from residuals
fit_resid <- lmFit(resid_expr, design_simple)
# contrast...
# ...continue with contrasts as above
