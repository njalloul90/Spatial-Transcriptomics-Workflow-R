# 01-qc-preprocess.R
# Run QC & preprocessing for multiple spatial samples and merge

source("scripts/00-setup.R")
source("scripts/utils.R")

# === user edits: list your samples and paths here ===
# Example:
# samples <- list(
#   "healthy_rep1" = "data/raw/healthy_rep1",
#   "disease_rep1" = "data/raw/disease_rep1",
#   ...
# )
samples <- list.files("data/raw", full.names = TRUE)
sample_names <- basename(samples)

seurat_list <- list()
for(i in seq_along(samples)){
  sd <- samples[i]
  name <- sample_names[i]
  # adapt depending on layout: Load10X_Spatial or Read10X
  s <- Load10X_Spatial(data.dir = sd)
  s$sample_id <- name
  # basic QC metrics
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  s[["nFeature_RNA"]] <- nFeature_RNA(s)
  s[["nCount_RNA"]] <- nCount_RNA(s)
  # filter: tune thresholds based on dataset
  s <- subset(s, subset = nFeature_RNA > 200 & percent.mt < 20)
  # Normalize & find variable features
  s <- SCTransform(s, assay="Spatial", verbose = FALSE) # or SCTransform on RNA assay
  s <- RunPCA(s, verbose = FALSE)
  s <- FindNeighbors(s, reduction = "pca", dims = 1:30)
  s <- FindClusters(s, resolution = 0.4)
  s <- RunUMAP(s, reduction = "pca", dims = 1:30)
  seurat_list[[name]] <- s
}

# optional: integrate samples (anchors) if you want a joint embedding
# Here we create a merged Seurat object for downstream joint clustering/annotation
# convert SCTransform outputs into a common assay if needed
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)
seurat_merged <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# run dimensionality reduction & clustering on merged
seurat_merged <- RunPCA(seurat_merged)
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:30)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.5)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:30)

saveRDS(seurat_merged, file = file.path(proc_dir, "seurat_merged_qc.rds"))
message("Saved preprocessed Seurat object to: ", file.path(proc_dir, "seurat_merged_qc.rds"))
