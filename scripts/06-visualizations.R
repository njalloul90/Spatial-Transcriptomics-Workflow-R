# 06-visualizations.R
source("scripts/00-setup.R")
seurat_annot <- readRDS(file.path(proc_dir, "seurat_merged_annotated.rds"))

# UMAP with cell types
p1 <- DimPlot(seurat_annot, reduction = "umap", group.by = "cell_type", label = TRUE) + ggtitle("UMAP: cell types")
ggsave(file.path(fig_dir, "umap_celltypes.png"), p1, width = 8, height = 6)

# spatial plots: example feature and cell type distribution
# pick top gene from earlier top markers file
top_markers <- read.csv(file.path(tables_dir, "cluster_top_markers.csv"))
example_gene <- top_markers$gene[1]
p2 <- SpatialFeaturePlot(seurat_annot, features = example_gene) + ggtitle(example_gene)
ggsave(file.path(fig_dir, paste0("spatial_", example_gene, ".png")), p2, width = 8, height = 6)

# plot cell type proportions across samples
prop_table <- read.csv(file.path(tables_dir, "celltype_proportions_by_sample.csv"))
p3 <- prop_table %>% ggplot(aes(x = sample_id, y = prop, fill = cell_type)) +
  geom_col(position = "fill") + theme_minimal() + coord_flip() +
  ggtitle("Cell type proportions by sample")
ggsave(file.path(fig_dir, "celltype_proportions_by_sample.png"), p3, width = 10, height = 6)
