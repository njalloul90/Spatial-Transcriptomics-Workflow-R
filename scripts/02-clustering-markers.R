# 02-clustering-markers.R
# Find cluster markers for the integrated object

source("scripts/00-setup.R")
seurat_merged <- readRDS(file.path(proc_dir, "seurat_merged_qc.rds"))

# use active idents (clusters)
Idents(seurat_merged) <- seurat_merged$seurat_clusters

# find markers per cluster (Wilcoxon)
markers <- FindAllMarkers(seurat_merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers, file.path(tables_dir, "cluster_markers_all.csv"), row.names = FALSE)

# Top markers per cluster
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top_markers, file.path(tables_dir, "cluster_top_markers.csv"), row.names = FALSE)

# Save a figure: spatial feature plots for top marker of each cluster (first marker)
pdf(file.path(fig_dir, "cluster_top_marker_spatial.pdf"), width = 9, height = 6)
for(k in unique(top_markers$cluster)){
  mk <- top_markers %>% filter(cluster == k) %>% slice(1) %>% pull(gene)
  p <- SpatialFeaturePlot(seurat_merged, features = mk, images = NULL) + ggtitle(paste0("Cluster ", k, " marker: ", mk))
  print(p)
}
dev.off()
