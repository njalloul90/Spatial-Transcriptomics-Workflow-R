# 03-celltype-annotation.R
# Two approaches: SingleR reference-based and manual marker-based annotation

source("scripts/00-setup.R")
seurat_merged <- readRDS(file.path(proc_dir, "seurat_merged_qc.rds"))

# convert to SingleCellExperiment for SingleR
sce <- as.SingleCellExperiment(seurat_merged, assay = "SCT")

# load a reference (e.g., HumanPrimaryCellAtlasData, MonacoImmuneData, or a custom reference)
ref <- celldex::HumanPrimaryCellAtlasData()  # choose as appropriate for organism/dataset

# run SingleR
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
# attach predictions back to Seurat metadata
seurat_merged$SingleR_label <- pred$labels

# Optionally smooth or propagate labels across spatial neighbors / clusters
seurat_merged$cluster_label <- paste0("cl", seurat_merged$seurat_clusters)
# Summarize per cluster to get cluster-level labels (majority vote)
cluster_labels <- data.frame(cluster = seurat_merged$seurat_clusters, SingleR = seurat_merged$SingleR_label)
cl_map <- cluster_labels %>% group_by(cluster) %>% count(SingleR) %>% top_n(1, n) %>% ungroup()
cl_map <- cl_map %>% select(cluster, SingleR)
# create a cluster->label map
cluster_to_label <- setNames(as.character(cl_map$SingleR), cl_map$cluster)
seurat_merged$cell_type <- as.character(cluster_to_label[as.character(seurat_merged$seurat_clusters)])
# fallback to SingleR label if cluster assignment NA
seurat_merged$cell_type[is.na(seurat_merged$cell_type)] <- seurat_merged$SingleR_label[is.na(seurat_merged$cell_type)]

# Save annotated Seurat object
saveRDS(seurat_merged, file = file.path(proc_dir, "seurat_merged_annotated.rds"))
write.csv(table(seurat_merged$cell_type, seurat_merged$sample_id), file.path(tables_dir, "celltype_by_sample_counts.csv"))
