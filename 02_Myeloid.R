## ============================================================
## AAA + LiCl
## Myeloid reclustering:
## Macrophage + Monocyte + Neutrophil
## ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

group_colors <- c(
  "50"  = "#FFD966B2",
  "100" = "#6FA8DCB2",
  "150" = "#8E7CC3B2",
  "200" = "#E06666B2"
)

set.seed(1234)

outdir <- path.expand("/Users/zhangshuli/Documents/R-training/Li/scRNA_li_AAA/03_myeloid")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("OUTPUT DIR:\n", outdir, "\n\n")

myeloid_celltypes <- c("Macrophage", "Monocyte", "Neutrophil")
myeloid <- subset(sce, subset = celltype %in% myeloid_celltypes)
myeloid$celltype_major <- myeloid$celltype


DefaultAssay(myeloid) <- "RNA"
## Seurat v5 objects may contain multiple RNA layers after merge
if (inherits(myeloid[["RNA"]], "Assay5")) {
  myeloid <- JoinLayers( myeloid, assay = "RNA")
}

myeloid <- NormalizeData(
  object = myeloid, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)

myeloid <- FindVariableFeatures(
  object = myeloid, selection.method = "vst", nfeatures = 3000, verbose = TRUE)

myeloid <- ScaleData(object = myeloid, features = VariableFeatures(myeloid), verbose = TRUE)
myeloid <- RunPCA(object = myeloid, features = VariableFeatures(myeloid), npcs = 50, verbose = TRUE)

p_elbow <- ElbowPlot(myeloid, ndims = 50) +
  theme_classic( base_size = 12, base_family = "sans")

ggsave(
  filename = file.path(outdir, "01_myeloid_PCA_elbow.pdf"),
  plot = p_elbow, width = 5, height = 4)

myeloid <- FindNeighbors(object = myeloid, reduction = "pca", dims = 1:30, verbose = TRUE)
myeloid <- FindNeighbors(myeloid, reduction = "harmony", dims = 1:30)
myeloid <- RunUMAP(myeloid, reduction = "harmony", dims = 1:30)

resolution_set <-seq(0.2, 1.0, 0.1)

for (res in resolution_set) {
  myeloid <- FindClusters(
    object = myeloid, resolution = res,
    algorithm = 1, random.seed = 1234, verbose = TRUE)
}

colnames(myeloid@meta.data)

clustree(myeloid@meta.data, prefix = "RNA_snn_res.")
ggsave(filename = "./clustree.pdf", width = 12, height = 10)

p1 <- DimPlot(myeloid, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T, repel = F, shuffle = T)
p2 <- DimPlot(myeloid, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T, repel = F, shuffle = T)
p3 <- DimPlot(myeloid, reduction = "umap", group.by = "RNA_snn_res.0.4", label = T, repel = F, shuffle = T)
p4 <- DimPlot(myeloid, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T, repel = F, shuffle = T)
p5 <- DimPlot(myeloid, reduction = "umap", group.by = "RNA_snn_res.0.6", label = T, repel = F, shuffle = T)
p6 <- DimPlot(myeloid, reduction = "umap", group.by = "RNA_snn_res.0.7", label = T, repel = F, shuffle = T)
p7 <- DimPlot(myeloid, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T, repel = F, shuffle = T)
p8 <- DimPlot(myeloid, reduction = "umap", group.by = "RNA_snn_res.0.9", label = T, repel = F, shuffle = T)
p9 <- DimPlot(myeloid, reduction = "umap", group.by = "RNA_snn_res.1", label = T, repel = F, shuffle = T)

p <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3)
ggsave(p, filename = "./RNA_snn_res.0.2-0.8_0930.pdf", width = 22, height = 20) # 分辨率为0.4,"RNA_snn_res.0.4"

Idents(myeloid) <- myeloid$RNA_snn_res.0.4

myeloid <- RunUMAP(
  object = myeloid,
  reduction = "pca",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.30,
  seed.use = 1234,
  verbose = TRUE
)

major_colors <- c(
  "Macrophage"    = "#BC3C29", 
  "Monocyte"      = "#38618B",
  "Neutrophil"    = "#5A7A7C")

NEJM_11 <- colorRampPalette(ggsci::pal_nejm()(8))(11)
show_col(NEJM_11)

p_major <- DimPlot(
  object = myeloid,
  reduction = "umap",
  group.by = "celltype_major",
  cols = major_colors,
  pt.size = 1.8,
  shuffle = TRUE,
  raster = TRUE) +
  coord_fixed() +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")

ggsave("02_myeloid_UMAP_major_celltype.pdf", p_major, width = 4.5, height = 5)

my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', 
                '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', 
                '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
                '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', 
                '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', 
                '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

p_cluster <- DimPlot(
  object = myeloid, reduction = "umap",
  group.by = "RNA_snn_res.0.4",
  label = TRUE, repel = TRUE,
  pt.size = 1.5, cols = my36colors,
  shuffle = TRUE, raster = FALSE) +
  coord_fixed() +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    axis.title = element_text(size = 12)) +
  labs(x = "UMAP 1", y = "UMAP 2")

ggsave("./03_myeloid_UMAP_clusters_res0.4.pdf",
  p_cluster, width = 4, height = 4)

DimPlot(
  object = myeloid, reduction = "umap",
  group.by = "RNA_snn_res.0.4",
  split.by = "sample",
  label = TRUE, label.size = 4.2, label.box = FALSE,
  repel = TRUE,
  pt.size = 1.5, 
  cols = my36colors,
  shuffle = TRUE, raster = FALSE) +
  coord_fixed() +
  labs(
    title = NULL,
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    ## axis
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(
      size = 12, colour = "black",
      margin = margin(t = 6)),
    axis.title.y = element_text(
      size = 12, colour = "black",
      margin = margin(r = 6)),
    axis.line = element_line(
      linewidth = 0.5, colour = "black"),
    ## facet title
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 12,
      face = "bold",
      colour = "black",
      margin = margin(
        t = 5,
        b = 7
      )
    ),
    
    ## legend
    legend.title = element_blank(),
    legend.text = element_text(size = 9, colour = "black"),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.key.width = grid::unit(0.35, "cm"),
    ## overall spacing
    panel.spacing.x = grid::unit(0.25, "cm"),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )
ggsave("./myeloid_UMAP_clusters_split.pdf", width = 11, height = 5)

Idents(myeloid) <- myeloid$RNA_snn_res.0.4

markers <- FindAllMarkers(
  object = myeloid,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, "./myeloid_cluster_markers_res0.4.csv", row.names = FALSE)

genes_to_check <- c(
  "Ly6c2", "Ccr2", "Plac8", "Vcan",    # Classical monocyte
  "Cd163", "Selenop", "Mrc1", "Folr2", # Resident macrophage
  "Spp1", "Mmp12", "Lgals3", "Arg1",   # Remodeling macrophage
  "Ctsk", "Atp6v0d2", "Oscar", "Acp5", # Osteoclast-like Mac
  "Camp", "Ngp", "Ltf",                # Immature neutrophil
  "Chil3", "S100a8", "S100a9", "Il1b", # Inflammatory Neu
  "Cd101", "Hcar2", "Cxcr4", "Retnlg", # Mature/CXCR4hi Neu
  "Gdf15", "Ccl8", "Ms4a7", # C0 monocyte-derived / transitional macrophage
  "Xcr1", "Clec9a",         # cDC1
  "Gp9", "Gp6", "Itga2b",    # Platelet aggregate
  "Hba-a1", "Hbb-bs", "Alas2",  # Hb-high QC
  "Ccr7", "Cd209a", "Il4i1"  # Migratory DC
  # "Col1a1", "Col1a2", "Col3a1" 
)
  
genes_to_check <- rev(unique(genes_to_check))

DotPlot(myeloid, features = genes_to_check, assay='RNA',
        cols = c("#F7FBFF", "#08306B")
        # cols = c("#F7FBFF", "#08519C")
        ) + coord_flip()

p_dot <- DotPlot(
  object = myeloid,
  group.by = "myeloid_subtype_short",
  features = genes_to_check,
  assay = "RNA",
  cols = c("#F7FBFF", "#08306B"),
  dot.scale = 6, scale = TRUE) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      size = 10, colour = "black",
      angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(
      size = 11, colour = "black"),  # , face = "italic"
    axis.ticks = element_blank(),
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 14, face = "bold", colour = "black"),
    legend.text = element_text(size = 14, colour = "black"),
    legend.key.height = grid::unit(0.45, "cm"),
    legend.key.width = grid::unit(0.45, "cm"),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 8 ) ) +
  guides(size = guide_legend(title = "Percent expressed", order = 1),
    colour = guide_colourbar(
      title = "Average expression",
      order = 2,
      barheight = grid::unit(3.5, "cm"),
      barwidth = grid::unit(0.45, "cm")
    )
  )
p_dot
ggsave("./myeloid_subtype_marker_DotPlot.pdf",
  plot = p_dot, width = 6, height = 7.5)


marker_modules <- list(
  "Classical monocyte" = c("Ly6c2", "Ccr2", "Plac8", "Vcan"),
  "Transitional macrophage" = c("Gdf15", "Ccl8", "Ms4a7", "Lpl"),
  "Resident macrophage" = c("Cd163", "Selenop", "Mrc1", "Folr2"),
  "Remodeling macrophage" = c("Spp1", "Mmp12", "Lgals3", "Arg1"),
  "Osteoclast-like macrophage" = c("Ctsk", "Atp6v0d2", "Oscar", "Acp5"),
  "Immature neutrophil" = c("Camp", "Ngp", "Ltf"),
  "Inflammatory neutrophil" = c("Chil3", "S100a8", "S100a9", "Il1b"),
  "Mature/CXCR4hi neutrophil" = c("Cd101", "Hcar2", "Cxcr4", "Retnlg"),
  "Migratory DC" = c("Ccr7", "Cd209a", "H2-DMb2", "Il4i1"),
  "cDC1" = c("Xcr1", "Clec9a", "Gcsam", "Mycl"),
  "Fibrocyte-like" = c("Ptprc", "Lyz2", "Col1a1", "Col1a2", "Col3a1"),
  "Hb-high" = c("Hba-a1", "Hbb-bs", "Alas2"),
  "Platelet-associated" = c("Gp9", "Gp6", "Itga2b", "Clec1b", "Treml1")
)

library(scRNAtoolVis)
pdf("./UMAP_split_myeloid_1.pdf", width = 13, height = 6)
featurePlot(object = myeloid,
            genes = c("Ly6c2", "Gdf15", "Cd163", "Spp1", "Ctsk", "Camp", "Chil3", "Cxcr4"),
            # genes =  c(  "Ccr2",     # Classical Mono
            #              "Ccl8",     # Transitional Mac
            #              "Selenop",  # Resident Mac
            #              "Mmp12",    # Remodeling Mac
            #              "Acp5",     # Osteoclast-like Mac
            #              "Ngp",      # Immature Neu
            #              "S100a8",   # Inflammatory Neu
            #              "Cd101"     # Mature/CXCR4hi Neu
            # ),
            quantile.val = 1,
            color = c("#F2F2F2", "#DCE6F2", "#AFC4DD", "#627FA8", "#111455"),
            nrow = 2, ncol = 4, add.corArrow = T)
dev.off()

library(scCustomize)
gene <- c("Ccr2",     # Classical Mono
          "Ccl8",     # Transitional Mac
          "Selenop",  # Resident Mac
          "Mmp12",    # Remodeling Mac
          "Acp5",     # Osteoclast-like Mac
          "Ngp",      # Immature Neu
          "S100a8",   # Inflammatory Neu
          "Cd101"     # Mature/CXCR4hi Neu
)
i=1
plots=list()
for (i in 1:length(gene)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = myeloid, reduction = "umap",
                                  colors_use = colorRampPalette(c("#FFF7F3", "#F768A1", "#49006A"))(50),
                                  features = gene[i]) + NoAxes() 
}
library(patchwork)
p <- wrap_plots(plots, ncol = 4)
p


marker_modules <- lapply(marker_modules, function(x) intersect(x, rownames(myeloid)))
marker_modules <- marker_modules[lengths(marker_modules) > 0]
genes_heatmap <- unique(unlist(marker_modules, use.names = FALSE))

avg_exp <- AverageExpression(
  object = myeloid, assays = "RNA",
  features = genes_heatmap,
  group.by = "myeloid_subtype_short",
  layer = "data", verbose = FALSE)$RNA

subtype_order_heatmap <- c(
  "Classical Mono",
  "Transitional Mac",
  "Resident Mac",
  "Remodeling Mac",
  "Osteoclast-like Mac",
  "Fibrocyte-like",
  "Immature Neu",
  "Inflammatory Neu",
  "Mature/CXCR4hi Neu",
  "Migratory DC",
  "cDC1",
  "Hb-high",
  "Platelet-associated"
)

subtype_order_heatmap <- intersect(subtype_order_heatmap, colnames(avg_exp))
avg_exp <- avg_exp[genes_heatmap, subtype_order_heatmap, drop = FALSE]

avg_z <- t(scale(t(avg_exp)))
avg_z[avg_z >  2.5] <-  2.5
avg_z[avg_z < -2.5] <- -2.5

gene_module <- rep(names(marker_modules), lengths(marker_modules))
names(gene_module) <- unlist(marker_modules, use.names = FALSE)
gene_module <- gene_module[rownames(avg_z)]
gene_module <- factor(gene_module, levels = names(marker_modules))

module_to_subtype <- c(
  "Classical monocyte"          = "Classical Mono",
  "Transitional macrophage"     = "Transitional Mac",
  "Resident macrophage"         = "Resident Mac",
  "Remodeling macrophage"       = "Remodeling Mac",
  "Osteoclast-like macrophage"  = "Osteoclast-like Mac",
  "Immature neutrophil"         = "Immature Neu",
  "Inflammatory neutrophil"     = "Inflammatory Neu",
  "Mature/CXCR4hi neutrophil"   = "Mature/CXCR4hi Neu",
  "Migratory DC"                = "Migratory DC",
  "cDC1"                        = "cDC1",
  "Fibrocyte-like"              = "Fibrocyte-like",
  "Hb-high"                     = "Hb-high",
  "Platelet-associated"         = "Platelet-associated"
)

row_title_colors <- subtype_colors[module_to_subtype[levels(gene_module)]]
row_name_colors <- subtype_colors[module_to_subtype[as.character(gene_module)]]

subtype_lineage <- c(
  "Classical Mono"        = "Monocyte",
  
  "Transitional Mac"      = "Macrophage",
  "Resident Mac"          = "Macrophage",
  "Remodeling Mac"        = "Macrophage",
  "Osteoclast-like Mac"   = "Macrophage",
  
  "Fibrocyte-like"        = "Fibrocyte-like",
  
  "Immature Neu"          = "Neutrophil",
  "Inflammatory Neu"      = "Neutrophil",
  "Mature/CXCR4hi Neu"    = "Neutrophil",
  
  "Migratory DC"          = "Dendritic cell",
  "cDC1"                  = "Dendritic cell",
  
  "Hb-high"               = "QC",
  "Platelet-associated"   = "QC"
)

lineage_colors <- c(
  "Monocyte"       = "#D39450",
  "Macrophage"     = "#5B78B8",
  "Fibrocyte-like" = "#8C6BB1",
  "Neutrophil"     = "#C56860",
  "Dendritic cell" = "#4E927A",
  "QC"             = "#9A9A9A"
)

ha_top <- HeatmapAnnotation(
  Lineage = subtype_lineage[colnames(avg_z)],
  col = list(Lineage = lineage_colors),
  annotation_name_gp = gpar(fontsize = 12, fontfamily = "sans", fontface = "bold"),
  simple_anno_size = unit(0.35, "cm")
)

column_name_colors <- subtype_colors[colnames(avg_z)]

col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("#F7FBFF", "#DCE6F2", "#111455"))

ht_marker <- Heatmap(
  avg_z, name = "Expression\nZ-score",
  col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  ## 保留分组和组间空隙
  row_split = gene_module,
  ## 但是不显示最左侧 cell subtype/module 名称
  row_title = NULL,
  top_annotation = ha_top,
  border = FALSE,
  rect_gp = gpar(col = "white", lwd = 0.8),
  ## Marker gene names 按对应 subtype 着色
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontfamily = "sans", col = row_name_colors),
  ## Bottom subtype names 按 subtype 着色
  column_names_gp = gpar(
    fontsize = 12, fontfamily = "sans",
    col = column_name_colors, fontface = "bold"),
  column_names_rot = 45,
  ## 仍然保留不同 marker module 之间的间距
  row_gap = unit(1.8, "mm"),
  column_gap = unit(1, "mm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(3, "cm")
  )
)

pdf("./myeloid_subtype_marker_heatmap_colored.pdf",
    width = 7, height = 11, family = "sans")
draw(
  ht_marker,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)
dev.off()

table(myeloid$RNA_snn_res.0.4, myeloid$celltype_major)
tab_origin <- table(myeloid$RNA_snn_res.0.4, myeloid$celltype_major)
round(prop.table(tab_origin, margin = 1) * 100, 2)

cluster_col <- "RNA_snn_res.0.4"
cluster_id <- as.character(myeloid@meta.data[[cluster_col]])

annotation_map <- c(
  "0"  = "Gdf15+Ccl8+ transitional macrophage",
  "1"  = "CXCR4hi mature neutrophil",
  "2"  = "Cd163+Selenop+ resident-like macrophage",
  "3"  = "CCR2+Ly6c2+ classical monocyte",
  "4"  = "Spp1+Mmp12+ remodeling macrophage",
  "5"  = "Chil3+ inflammatory neutrophil",
  "6"  = "CCR7+Cd209a+ dendritic cell",
  "7"  = "COL1A1+ fibrocyte-like myeloid",
  "8"  = "Hb-high erythroid-associated cells",
  "9"  = "Xcr1+Clec9a+ cDC1",
  "10" = "Camp+Ngp+ immature neutrophil",
  "11" = "Platelet-associated cells",
  "12" = "Ctsk+Atp6v0d2+ osteoclast-like macrophage"
)

myeloid$myeloid_subtype <- unname(annotation_map[cluster_id])
table(myeloid$myeloid_subtype)

lineage_map <- c(
  "0"  = "Macrophage",
  "1"  = "Neutrophil",
  "2"  = "Macrophage",
  "3"  = "Monocyte",
  "4"  = "Macrophage",
  "5"  = "Neutrophil",
  "6"  = "Dendritic_cell",
  "7"  = "Fibrocyte-like_myeloid",
  "8"  = "QC_Hb-high",
  "9"  = "Dendritic_cell",
  "10" = "Neutrophil",
  "11" = "QC_Platelet",
  "12" = "Macrophage"
)

myeloid$myeloid_lineage <- unname(lineage_map[cluster_id])

core_clusters <- c("0", "1", "2", "3", "4", "5", "7", "10", "12")
myeloid$core_myeloid <- cluster_id %in% core_clusters
exclude_clusters <- c("8", "11")
myeloid$exclude_from_downstream <- cluster_id %in% exclude_clusters
table(Cluster = cluster_id, Subtype = myeloid$myeloid_subtype)
table(myeloid$myeloid_lineage)
table(myeloid$core_myeloid)
table(myeloid$exclude_from_downstream)


cluster_col <- "RNA_snn_res.0.4"
cluster_id <- as.character(myeloid@meta.data[[cluster_col]])
short_annotation_map <- c(
  "0"  = "Transitional Mac",
  "1"  = "Mature/CXCR4hi Neu",
  "2"  = "Resident Mac",
  "3"  = "Classical Mono",
  "4"  = "Remodeling Mac",
  "5"  = "Inflammatory Neu",
  "6"  = "Migratory DC",
  "7"  = "Fibrocyte-like",
  "8"  = "Hb-high",
  "9"  = "cDC1",
  "10" = "Immature Neu",
  "11" = "Platelet-associated",
  "12" = "Osteoclast-like Mac"
)

myeloid$myeloid_subtype_short <- unname(short_annotation_map[cluster_id])
short_order <- c(
  "Classical Mono",
  "Transitional Mac",
  "Resident Mac",
  "Remodeling Mac",
  "Osteoclast-like Mac",
  "Fibrocyte-like",
  "Immature Neu",
  "Inflammatory Neu",
  "Mature/CXCR4hi Neu",
  "Migratory DC",
  "cDC1",
  "Hb-high",
  "Platelet-associated"
)

myeloid$myeloid_subtype_short <- factor(myeloid$myeloid_subtype_short,levels = short_order)

subtype_colors <- c(
  ## Monocyte
  "Classical Mono" = "#D98C4A",
  ## Macrophage
  "Transitional Mac" = "#AFD8E2", "Resident Mac" = "#7ED0F0", 
  "Remodeling Mac" = "#5A61A8",   "Osteoclast-like Mac" = "#111455",
  ## Fibrocyte-like
  "Fibrocyte-like" = "#8C6BB1",
  ## Neutrophil
  "Immature Neu" = "#E7A9A1",
  "Inflammatory Neu" = "#D4675F",
  "Mature/CXCR4hi Neu" = "#A83E3E",
  ## DC
  "Migratory DC" = "#65A98F", "cDC1" = "#2F7F65",
  ## QC / associated
  "Hb-high" = "#BEBEBE", "Platelet-associated" = "#8C8C8C"
)

p_myeloid_annotation <- DimPlot(
  object = myeloid,
  reduction = "umap",
  group.by = "myeloid_subtype_short",
  cols = subtype_colors,
  label = FALSE,
  repel = TRUE,
  label.size = 4,
  label.box = FALSE,
  pt.size = 0.45,
  shuffle = TRUE,
  raster = FALSE
) +
  coord_fixed() +
  labs(x = "UMAP 1",y = "UMAP 2", title = NULL) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 11, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10, colour = "black"),
    legend.key.height = grid::unit(0.4, "cm"),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )

p_myeloid_annotation
ggsave("./myeloid_sub_celltype.pdf", p_myeloid_annotation, width = 5, height = 5)

DimPlot(
  object = myeloid, reduction = "umap",
  group.by = "myeloid_subtype_short",
  split.by = "sample",
  label = FALSE, 
  # label.size = 4.2, label.box = FALSE,
  repel = TRUE,
  pt.size = 0.5, 
  cols = subtype_colors,
  shuffle = TRUE, raster = FALSE) +
  coord_fixed() +
  labs(
    title = NULL,
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    ## axis
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(
      size = 12, colour = "black",
      margin = margin(t = 6)),
    axis.title.y = element_text(
      size = 12, colour = "black",
      margin = margin(r = 6)),
    axis.line = element_line(
      linewidth = 0.5, colour = "black"),
    ## facet title
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 12,
      face = "bold",
      colour = "black",
      margin = margin(
        t = 5,
        b = 7
      )
    ),
    
    ## legend
    legend.title = element_blank(),
    legend.text = element_text(size = 9, colour = "black"),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.key.width = grid::unit(0.35, "cm"),
    ## overall spacing
    panel.spacing.x = grid::unit(0.25, "cm"),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )
ggsave("./myeloid_UMAP_subcelltype_split.pdf", width = 11, height = 5)


core_subtypes <- c(
  "Classical Mono",
  "Transitional Mac",
  "Resident Mac",
  "Remodeling Mac",
  "Osteoclast-like Mac",
  "Fibrocyte-like",
  "Immature Neu",
  "Inflammatory Neu",
  "Mature/CXCR4hi Neu"
)

macrophage_subtypes <- c(
  "Transitional Mac",
  "Resident Mac",
  "Remodeling Mac",
  "Osteoclast-like Mac"
)

neutrophil_subtypes <- c(
  "Immature Neu",
  "Inflammatory Neu",
  "Mature/CXCR4hi Neu"
)

meta_core <- myeloid@meta.data %>%
  mutate(
    sample = as.character(sample),
    subtype = as.character(myeloid_subtype_short)
  ) %>% filter(subtype %in% core_subtypes)

core_comp <- meta_core %>%
  count(sample, subtype, name = "n_cells") %>%
  complete(
    sample = c("50", "100", "150", "200"),
    subtype = core_subtypes,
    fill = list(n_cells = 0)) %>%
  group_by(sample) %>%
  mutate(
    total_selected_myeloid = sum(n_cells),
    percentage = n_cells / total_selected_myeloid * 100
  ) %>% ungroup()
core_comp$sample <- factor(as.character(core_comp$sample), levels = c("50", "100", "150", "200"))
write.csv(core_comp, "./01_core_myeloid_subtype_composition.csv", row.names = FALSE)

mac_comp <- meta_core %>%
  filter(subtype %in% macrophage_subtypes) %>%
  count(sample, subtype, name = "n_cells") %>%
  group_by(sample) %>%
  mutate(
    total_macrophage = sum(n_cells),
    pct_within_macrophage = n_cells / total_macrophage * 100
  ) %>% ungroup()
write.csv(mac_comp, "./02_macrophage_state_composition.csv", row.names = FALSE)

neu_comp <- meta_core %>%
  filter(subtype %in% neutrophil_subtypes) %>%
  count(sample, subtype, name = "n_cells") %>%
  group_by(sample) %>%
  mutate(
    total_neutrophil = sum(n_cells),
    pct_within_neutrophil = n_cells / total_neutrophil * 100
  ) %>% ungroup()
write.csv(neu_comp, "03_neutrophil_state_composition.csv", row.names = FALSE)

core_comp$sample <- factor(as.character(core_comp$sample), levels = c("50", "100", "150", "200"))

subtype_order <- c(
  "Classical Mono",
  "Transitional Mac",
  "Resident Mac",
  "Remodeling Mac",
  "Osteoclast-like Mac",
  "Fibrocyte-like",
  "Immature Neu",
  "Inflammatory Neu",
  "Mature/CXCR4hi Neu"
)

core_comp$subtype <- factor(as.character(core_comp$subtype), levels = rev(subtype_order))

lineage_map <- c(
  "Classical Mono"        = "Monocyte",
  "Transitional Mac"      = "Macrophage",
  "Resident Mac"          = "Macrophage",
  "Remodeling Mac"        = "Macrophage",
  "Osteoclast-like Mac"   = "Macrophage",
  "Fibrocyte-like"        = "Fibrocyte-like",
  "Immature Neu"          = "Neutrophil",
  "Inflammatory Neu"      = "Neutrophil",
  "Mature/CXCR4hi Neu"    = "Neutrophil"
)

core_comp$lineage <- unname(
  lineage_map[as.character(core_comp$subtype)]
)

core_comp$lineage <- factor(
  core_comp$lineage, levels = c("Neutrophil", "Fibrocyte-like", "Macrophage", "Monocyte"))

core_comp <- core_comp %>% dplyr::mutate(percentage_label = sprintf("%.1f", percentage))

p_comp_heatmap <- ggplot(
  core_comp,
  aes(x = sample, y = subtype, fill = percentage)) +
  geom_tile(colour = "white", linewidth = 1.1, width = 0.96, height = 0.94) +
  # Percentage directly shown in each cell
  geom_text(
    aes(label = percentage_label, colour = percentage),
    size = 4.0, family = "sans", fontface = "plain") +
  ## Blue gradient
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#DCEAF4", "#9ECAE1", "#5B8CC0", "#344C8C", "#111455"),
    name = "Percentage (%)") +
  ## Automatically use white text in high-value cells
  scale_colour_gradient2(
    low = "black", mid = "black",
    high = "white",
    midpoint = max(core_comp$percentage, na.rm = TRUE) * 0.58,
    guide = "none") +
  facet_grid(lineage ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "LiCl dose (mg/kg)", y = NULL) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    ## Remove background grid
    panel.grid = element_blank(),
    ## Full outer frame
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.65),
    ## Dose labels
    axis.text.x = element_text(size = 12, colour = "black", face = "plain"),
    ## Cell subtype labels
    axis.text.y = element_text(size = 11, colour = "black", hjust = 1),
    axis.title.x = element_text(size = 13, colour = "black", margin = margin(t = 10)),
    axis.ticks = element_blank(),
    ## Lineage labels
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 0, size = 10, face = "bold", colour = "black", margin = margin(r = 8)),
    ## Legend
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = grid::unit(1, "cm"),
    legend.key.width = grid::unit(0.42,  "cm"),
    ## Spacing between lineages
    panel.spacing.y = grid::unit(0.10, "cm"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10 )
  )
p_comp_heatmap
ggsave("./myeloid_percentage_heatmap.pdf", p_comp_heatmap, width = 7, height = 5)

p_comp_dot <- ggplot(
  core_comp,
  aes(x = sample, y = subtype)) +
  geom_point(aes(size = percentage, colour = percentage), alpha = 0.95) +
  scale_size_continuous(range = c(1.5, 7), name = "Percentage (%)") +
  scale_colour_gradientn(
    colours = c("#DCEAF4", "#7ED0F0", "#5B78B8", "#292F66", "#111455"),
    name = "Percentage (%)") +
  facet_grid(lineage ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "LiCl dose (mg/kg)", y = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.title.x = element_text(size = 13, margin = margin(t = 8)),
    axis.ticks = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, size = 10, face = "bold", margin = margin(r = 8)),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.spacing.y = grid::unit(0.1, "cm")
  )

p_comp_dot
ggsave("./myeloid_percentage_dotplot.pdf", p_comp_dot, width = 7, height = 5)


core_colors <- subtype_colors[core_subtypes]
core_comp$subtype <- factor(core_comp$subtype, levels = core_subtypes)
p_core_stack <- ggplot(
  core_comp,
  aes(x = sample, y = percentage, fill = subtype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = core_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))) +
  labs(
    x = "LiCl dose (mg/kg)",
    y = "Percentage of myeloid cells (%)",
    fill = NULL) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8
    ),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 11, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    axis.line = element_blank(),
    #axis.line = element_line(colour = "black",  linewidth = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, colour = "black"),
    legend.key.height = grid::unit(0.4, "cm"),
    plot.margin = margin(8, 8, 8, 8)
  )
p_core_stack
ggsave("./myeloid_percentage.pdf", p_core_stack, width = 5, height = 4.5)

mac_comp$subtype <- factor(mac_comp$subtype, levels = macrophage_subtypes)
mac_colors <- subtype_colors[macrophage_subtypes]
p_mac_stack <- ggplot(
  mac_comp,
  aes(x = sample, y = pct_within_macrophage, fill = subtype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = mac_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))) +
  labs(x = "LiCl dose (mg/kg)", y = "Percentage within macrophages (%)", fill = NULL) +
  theme_classic( base_size = 12, base_family = "sans") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 11,  colour = "black"),
    axis.line = element_blank(),
    # axis.line = element_line(linewidth = 0.5, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10, colour = "black")
  )
p_mac_stack
ggsave("./macrophage_percentage.pdf", p_mac_stack, width = 5, height = 4.5)


neu_comp$subtype <- factor(neu_comp$subtype, levels = neutrophil_subtypes)
neu_colors <- subtype_colors[neutrophil_subtypes]
p_neu_stack <- ggplot(
  neu_comp,
  aes(x = sample, y = pct_within_neutrophil, fill = subtype)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = neu_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))) +
  labs(x = "LiCl dose (mg/kg)", y = "Percentage within neutrophils (%)", fill = NULL) +
  theme_classic( base_size = 12, base_family = "sans") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 11,  colour = "black"),
    axis.line = element_blank(),
    # axis.line = element_line(linewidth = 0.5, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10, colour = "black")
  )
p_neu_stack
ggsave("./neutrophils_percentage.pdf", p_neu_stack, width = 5, height = 4.5)

saveRDS(myeloid, "./myeloid.rds")



