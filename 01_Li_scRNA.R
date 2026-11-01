rm(list=ls())
options(stringsAsFactors = F) 

library(clusterProfiler)
library(clustree)
library(cowplot)
library(cols4all)
library(colorblindcheck)
library(dplyr)
library(DOSE)
library(enrichplot)
library(EnhancedVolcano)
library(fgsea)
library(ggplot2)
library(ggalluvial)
library(ggstatsplot)
library(ggsci)
library(gridExtra)
library(gtable)
library(GSVA)
library(harmony)
library(kableExtra)
library(msigdbr)
library(Matrix)
library(org.Hs.eg.db)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(Seurat)
library(shinyjs)
library(tidyverse)
library(tinyarray)

cluster_colors <- c(
  "Adipocyte"     = "#1B1919B2",
  "B_cell"        = "#008B45B2",
  "Endothelial"   = "#BB0021B2",
  "Fibroblast"    = "#A20056B2", 
  "Lymphatic_EC"  = "#808180B2",
  "Macrophage"    = "#3B4992B2", 
  "Monocyte"      = "#5F559BB2", 
  "Neutrophil"    = "#631879B2", 
  "Proliferating_cell" = "#EFC000B2", 
  "Smooth_muscle_cell" = "#EE0000B2", 
  "T_cell"        = "#008280B2"        
)

celltype_colors <- cluster_colors

celltype_colors <- c(
  "#3B4992B2", "#A20056B2", "#EE0000B2", "#5F559BB2", "#008280B2", 
  "#BB0021B2", "#808180B2", "#EFC000B2", "#631879B2", "#008B45B2", "#1B1919B2"
)

group_colors <- c(
  "50"  = "#FFD966B2",
  "100" = "#6FA8DCB2",
  "150" = "#8E7CC3B2",
  "200" = "#E06666B2"
)

sce$sample <- factor(sce$sample, levels = c("50", "100", "150", "200"))
table(sce$sample)   # 输出顺序即为指定顺序
sample_colors <- c(
  "50"  = "#FFD966B2",
  "100" = "#6FA8DCB2",
  "150" = "#8E7CC3B2",
  "200" = "#E06666B2"
)

################## step1:导入数据 ##################  
setwd("/Users/zhangshuli/Documents/R-training/Li/scRNA_li_AAA/")
getwd()

dir <- "./scRNA_Li_AAA"
samples <- c("50", "100", "150", "200")
required_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")

file_check <- lapply(samples, function(pro) {
  sample_dir <- file.path(dir, pro)
  
  data.frame(
    sample = pro,
    file = required_files,
    exists = file.exists(file.path(sample_dir, required_files))
  )
})

file_check <- bind_rows(file_check)

# 每个样本创建一个Seurat对象
sceList <- lapply(samples, function(pro) {
  print(paste0("正在读取样本：", pro))
  counts <- Read10X(data.dir = file.path(dir, pro))
  # 如果10X文件包含多种feature type，只保留Gene Expression
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      stop(paste0("样本 ", pro, " 未找到Gene Expression矩阵。"))
    }
  }
  sce <- CreateSeuratObject(
    counts = counts,
    project = pro,
    min.cells = 5,
    min.features = 300
  )
  sce$sample <- pro
  sce$group <- pro
  sce$dose <- as.numeric(pro)
  return(sce)
})
names(sceList) <- samples
lapply(sceList, dim)
data.frame(
  sample = samples,
  n_gene = sapply(sceList, nrow),
  n_cell = sapply(sceList, ncol)
)

sce <- merge(
  x = sceList[[1]],
  y = sceList[-1],
  add.cell.ids = samples,
  project = "Li_AAA"
)
# Seurat v5合并后可能保留多个counts layer
if (packageVersion("SeuratObject") >= package_version("5.0.0")) {
  sce <- JoinLayers(sce)
}
sce$group <- factor(sce$group, levels = samples)
sce$dose <- as.numeric(as.character(sce$group))

sce
table(sce$orig.ident)
table(sce$group)
saveRDS(sce, "./00_RAW/sce_raw.rds")

sce <- readRDS("./00_RAW/sce_raw.rds")
as.data.frame(sce@assays$RNA@layers$counts[1:10, 1:2])
head(sce@meta.data, 10)
tail(sce@meta.data, 10)
table(sce@meta.data$orig.ident) 
# 100   150   200    50 
# 16051 13321 16033 14961 
#分组

################## step2:质控 ##################  
mito_genes <- rownames(sce)[grep("^MT-|^mt-", rownames(sce))] 
mito_genes # 13线粒体基因
# 计算线粒体百分比
sce <- PercentageFeatureSet(sce, pattern = "^MT-|^mt-", col.name = "percent_mito")
# sce.all[["percent.mito"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
fivenum(sce@meta.data$percent_mito)
sce$orig.ident <- factor(sce$orig.ident, levels = c("50", "100", "150", "200"))
table(sce$orig.ident)
# 检查数据的可用性
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), 
        group.by = "orig.ident",
        ncol = 3, cols = brewer.pal(9,"Set1")[1:6])
ggsave(filename = "./00_RAW/RAW_Vlnplot_feature_count_mito.pdf", width = 20, height = 8)

plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent_mito", 
                        group.by = "orig.ident", cols = brewer.pal(9,"Set1")[1:6])
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                        group.by = "orig.ident", cols = brewer.pal(9,"Set1")[1:6])
plot1 + plot2
ggsave(filename = "./00_RAW/RAW_FeatureScatter_feature_count_mito.pdf", width = 11, height = 5)

# >25% of reads from mitochondria genes /Use
# removal of cells with <200 or >2,500 genes detected
# resulting in ∼8,000 cells.
sce1 <- subset(sce, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & 
                 nCount_RNA >= 500 & nCount_RNA <= 30000 &
                 percent_mito < 20 )
VlnPlot(sce1, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), 
        group.by = "orig.ident",
        ncol = 3, cols = brewer.pal(9,"Set1")[1:6])
ggsave(filename = "./01_QC/QC_Vlnplot_feature_count_mito.pdf", width = 20, height = 8)

plot1 <- FeatureScatter(sce1, feature1 = "nCount_RNA", feature2 = "percent_mito", 
                        group.by = "orig.ident", cols = brewer.pal(9,"Set1")[1:6])
plot2 <- FeatureScatter(sce1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                        group.by = "orig.ident", cols = brewer.pal(9,"Set1")[1:6])
plot1 + plot2
ggsave(filename = "./01_QC/QC_FeatureScatter_feature_count_mito.pdf", width = 11, height = 5)

saveRDS(sce1, "./01_QC/sce_QC.rds")
sce <- readRDS("./01_QC/sce_QC.rds")

################## step3:多样本整合，降维聚类分群harmony ################## 
# 归一化
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
# 寻找高变异度基因
sce <- FindVariableFeatures(sce)
# 尺度变换
sce <- ScaleData(sce)
# 主成分分析
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
# 不进行批次矫正
sce <- RunUMAP(sce, reduction = "pca", dims = 1:10, 
               reduction.name = "umap_naive", reduction.key = 'naiveUMAP_')
# 使用harmony进行批次矫正
# 默认reduction = "pca"
# group.by.vars参数输入的是批次信息(group.by.vars = "group")
# reduction.save 是结果保存的名称
sce <- RunHarmony(sce, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:10, 
               reduction.name = "umap_harmony", reduction.key = 'harmonyUMAP_')
p1 <- DimPlot(sce, reduction = "umap_naive", label=T, cols = my36colors) # cols = brewer.pal(12,"Paired")[1:8]
p1_1 <- DimPlot(sce, reduction = "umap_harmony", label=T, cols = my36colors)
p1+p1_1
ggsave(p1+p1_1, filename="./01_QC/Dimplot_harmony_afterQC.pdf", width=12, height = 5)

################## step4: 降维聚类分群 ################## 
setwd("./01_QC/")
sce <- readRDS("./sce_QC.rds")
ElbowPlot(sce, ndims = 50, reduction = "harmony")
### 测试dims 
# 先默认分辨率设置为resolution=0.2
for (i in c(15, 20, 25, 30, 35, 40)) {
  sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.2)
  sce <- RunUMAP(sce, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(sce, reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p <- plot_15 + plot_20 + plot_25 + plot_30 + plot_35 + plot_40  + plot_layout(ncol = 3)
p
ggsave(p, filename = "./dim.15-40.pdf", width = 15, height = 8)
# 选择dims=15
### 测试resolution,设置不同的分辨率，观察分群效果 
# 设置一个随机种子，避免每次的图形是随机的
set.seed(123) 
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:15)
sce <- FindClusters(sce, resolution = seq(0.2, 0.8, 0.1))
### 多个分辨率的分群信息会保存在metadata中
clustree(sce@meta.data, prefix = "RNA_snn_res.")
ggsave(filename = "./clustree.pdf", width = 12, height = 10)

p1 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T, repel = F, shuffle = T)
p2 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T, repel = F, shuffle = T)
p3 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.4", label = T, repel = F, shuffle = T)
p4 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T, repel = F, shuffle = T)
p5 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.6", label = T, repel = F, shuffle = T)
p6 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.7", label = T, repel = F, shuffle = T)
p7 <- DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T, repel = F, shuffle = T)
p <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + plot_layout(ncol = 4)
ggsave(p, filename = "./RNA_snn_res.0.2-0.8_0930.pdf", width = 22, height = 10) # 分辨率为0.4,"RNA_snn_res.0.4"

# 查看特定分辨率的分群效果
colnames(sce@meta.data)
table(sce$RNA_snn_res.0.2)

# 删除一些RNA_snn_res
sce@meta.data <- sce@meta.data[ , -which(colnames(sce@meta.data) %in% c(
  "RNA_snn_res.0.3", "RNA_snn_res.0.4", "RNA_snn_res.0.5",
  "RNA_snn_res.0.6","RNA_snn_res.0.7", "RNA_snn_res.0.8"))]

# 最终降纬参数：dim=15, res=0.2
sce@reductions$umap <- NULL
sce@reductions$umap_test <- NULL

sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:15)
sce <- FindClusters(sce, resolution = 0.7)
table(sce@meta.data$seurat_clusters)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:15, 
                 reduction.name = 'umap', reduction.key = 'UMAP_')
sel.clust = "RNA_snn_res.0.7"
# sel.clust = "celltype"
sce <- SetIdent(sce, value = sel.clust)
table(sce@active.ident) 
ncluster <- length(unique(sce[[]]$RNA_snn_res.0.7))
mycol <- colorRampPalette(brewer.pal(12, "Paired"))(ncluster)
plot1 <- DimPlot(sce, reduction = "umap", cols = mycol, label = T, repel = F, 
                 pt.size = 2, shuffle = T)
ggsave(plot1, filename = "./umap_0830.pdf", width = 4.5, height = 4) 

################## step5:简单检查常见分群情况 ################## 
setwd("../02_Celltype/")

DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.3", pt.size = 2, cols = mycol, label = T) 
ggsave('umap_by_RNA_snn_res.0.3.pdf',width = 8,height = 6)

table(sce@meta.data$orig.ident)
DimPlot(sce, reduction = "umap", split.by = "orig.ident", cols = mycol, pt.size = 2, label = TRUE)
ggsave(filename = './UMAP_multi_samples_split.pdf', width = 16, height = 5)

# 查看细胞周期影响
FeaturePlot(sce,features = c("Pcna"), reduction = "umap", cols = brewer.pal(9, name = "RdPu"))
##s期基因PCNA，G2/M期基因MKI67均高表达，说明细胞均处于细胞周期
ggsave(file = "Cycling_Pcna_expression.pdf", width=6, height = 5)

# 鉴定marker基因
obj.combined <- sce
sample.markers <- FindAllMarkers(obj.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# 保存每个cluster top10的marker基因
top10_table <- unstack(top10, gene ~ cluster)
names(top10_table) <- gsub("X", "cluster", names(top10_table))
write.csv(file="top10_marker_genes.csv", top10_table, row.names=F)

genes_to_check <- c(
  'Myh11', 'Cnn1', 'Tagln', 'Acta2', 'Myl9',   
  # SMC, 'Tpm2', 'Smtn', 'Des', 'Myocd', 'Notch3',
  'Pecam1', 'Cdh5', 'Esam', 'Egfl7', 'Kdr', 
  # Endothelial,  'Vwf', 'Flt1', 'Ramp2', 'Emcn', 'Tek', 'Klf2',
  'Col1a1', 'Col1a2', 'Col3a1', 'Dcn', 
  # Fibroblast,  'Pdgfra', 'Lum', 'Col14a1', 'Pi16', 'Fbn1', 'Mfap5',
  'Prox1', 'Lyve1', 'Pdpn', 'Flt4', 
  # Lymphatic endothelial cells, 'Mmrn1', 'Ccl21a', 'Reln', 
  # Pericyte, 'Notch3', 'Des', 'Kcnj8', 'Rgs5', 'Pdgfrb', 'Cspg4', 'Abcc9', 'Mcam', 'Cald1',
  # 肌成纤维细胞, 'Postn', 'Cthrc1', 'Fn1', 'Col1a1', 'Acta2', 'Tagln', 'Timp1', 'Thbs1', 'Sparc',  
  # 神经鞘/施万细胞, 'Ncmap','Sox10', 'S100b', 'Plp1', 'Mpz', 'Mbp', 'Cnp',          
  # 间皮细胞,'Upk3b', 'Wt1', 'Krt19', 'Krt8', 'Krt18', 'Muc16', 'Msln',       
  'Adipoq', 'Plin1', 'Fabp4', 'Cidec', 'Lpl',     # 脂肪细胞, 'Pparg', 'Retn',
  'Ptprc', # CD45, 泛免疫细胞
  'C1qa' , 'C1qb', 'C1qc', 'Adgre1', 'Csf1r', 'Cd68', 
  # 巨噬细胞,'Fcgr1', 'Lgals3', 'Aif1', 'Lyz2', 
  'Lyz2', 'Ly6c2', 'Ccr2', 'Ctss', 
  # 单核细胞,'S100a8', 'S100a9', 'Fcgr3', 'Cx3cr1', 'Lst1', 
  'Csf3r', 'S100a8', 'S100a9',  
  # 中性粒细胞,'Ngp', 'Camp', 'Ly6g', 'Retnlg', 'Mmp8', 'Lcn2', 
  'H2-Ab1', 'Cd74', # 树突状细胞 DC,'Flt3', 'Itgax', 
  # cDC1, 'Xcr1', 'Clec9a', 'Batf3', 'Irf8', 'Cadm1', 
  # cDC2,'Mgl2', 'Clec10a', 'Cd209a', 'Sirpa', 'Irf4',  
  # 浆细胞样 DC,'Siglech', 'Bst2', 'Tcf4', 'Gzmb', 'Irf7', 
  'Cd3d', 'Cd3e', 'Cd3g', 'Trbc1',  
  # T 细胞, 'Lck', 'Lat', 'Trac',
  # CD4 T 细胞,'Foxp3', 'Il2ra', 'Ctla4', 'Malat1', 'Cd4', 'Il7r', 'Ltb', 
  # CD8 T 细胞,'Gzmb', 'Gzmk', 'Prf1', 'Cd8a', 'Cd8b1', 'Ccl5', 'Nkg7',   
  # NK 细胞,'Nkg7', 'Ncr1', 'Tyrobp', 'Klrd1', 'Prf1', 'Gzmb', 'Ccl5', 'Klrk1', 
  'Cd79a', 'Cd79b', 'Ms4a1', 'Cd37', 'Cd74', 
  # B 细胞, 'H2-Aa', 'Cd22', 'Cd19', 
  # 浆细胞, 'Derl3', 'Igha', 'Ighg', 'Sdc1', 'Xbp1', 'Jchain', 'Mzb1',
  # 'Kit', 'Cpa3', 'Mcpt4', 'Mcpt5', 'Tpsb2', 'Fcer1a', 'Ms4a2', # 肥大细胞
  # 'Hba-a1', 'Hba-a2', 'Hbb-bs', 'Hbb-bt', 'Alas2', 'Slc4a1', 'Gypa', # 红细胞/红系细胞
  # 血小板/巨核细胞,'Itga2b', 'Gp9', 'Gp1ba', 'Nrg1', 'Ppbp', 'Pf4',
  'Mki67', 'Top2a', 'Ube2c', 'Birc5'       # 增殖细胞	'Cenpf', 
  )    
# genes_to_check=str_to_upper(unique(genes_to_check)) #转换为大写
genes_to_check <- rev(unique(genes_to_check))
# th <- theme(axis.text.x = element_text(angle = 90, hjust=1))
p_all_markers <- DotPlot(sce, features = genes_to_check, assay='RNA',
                         cols = c("#F7FBFF", "#08519C")) + coord_flip() # + th 
p_all_markers
ggsave(p_all_markers, filename = "./markers to check_0830.pdf", width = 8, height = 13)

library(scRNAtoolVis)
cluster_order <- levels(Idents(sce))
cluster_colors_use <- cluster_colors[cluster_order]
pdf("./marker_heatmaps.pdf", width = 5.5, height = 8.5)
averageHeatmap(object = sce, markerGene = genes_to_check,
               htCol = c("#22496A", "#FDFDFD", "#8E348B"),
               # htCol = c("#22496A", "#FDFDFD", "#8E348B"),
               # htCol = c("#4E62AC", "white", "#9D0141"),
               cluster_rows = TRUE, 
               # cluster_columns = TRUE,
               annoCol = TRUE, myanCol = unname(cluster_colors_use))
dev.off()

################## step6:分群情况可视化 ################## 
sce
p_umap <- DimPlot(sce, reduction = "umap", group.by = "orig.ident",
                  pt.size = 1.5, cols = group_colors,
                  label = F, label.box = F)  
p_all_markers+p_umap
ggsave(p_umap, 'markers_umap_3.pdf',width = 15,height = 8)

DimPlot(sce, reduction = "umap",
        group.by = 'group',
        label = T) 


p_umap <- DimPlot(
  object = sce, reduction = "umap", group.by = "orig.ident",
  cols = group_colors,
  pt.size = 1.5,
  shuffle = TRUE,  raster = TRUE, label = FALSE) +
  coord_fixed() +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    axis.title = element_text(
      size = 12,
      face = "plain",
      colour = "black"
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, colour = "black"),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.45, "cm"),
    plot.title = element_text(
      size = 14, face = "bold",
      hjust = 0.5, colour = "black" ),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)) +
  labs(x = "UMAP 1", y = "UMAP 2")
p_umap

ggsave(filename = "UMAP_by_LiCl_dose.pdf",
  plot = p_umap,
  width = 6,
  height = 5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "UMAP_by_LiCl_dose.png",
  plot = p_umap,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)

ggsave('group_umap.pdf', width = 10, height = 9)

################## step7:细胞亚群命名 ################## 
# 0, 6, 11, 14: fibroblast
# 1: macrophage
# 2: monocyte
# 3: neutrophil
# 4,5: SMC
# 7 : endothelial
# 12: Lymphatic endothelial
# 8: T cell
# 9: Proliferating cell
# 10: B cell
# 13: adipocyte

celltype=data.frame(ClusterID=0:14, celltype= 0:14) 
#定义细胞亚群 
celltype[celltype$ClusterID %in% c(0, 6, 11, 14),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(7),2]='Endothelial'
celltype[celltype$ClusterID %in% c(12),2]='Lymphatic_EC'
celltype[celltype$ClusterID %in% c(4, 5),2]='Smooth_muscle_cell' 
celltype[celltype$ClusterID %in% c(2),2]='Monocyte' 
celltype[celltype$ClusterID %in% c(1),2]='Macrophage'
celltype[celltype$ClusterID %in% c(3),2]='Neutrophil'
celltype[celltype$ClusterID %in% c(8),2]='T_cell' 
celltype[celltype$ClusterID %in% c(9),2]='Proliferating_cell'
celltype[celltype$ClusterID %in% c(10),2]='B_cell'  
celltype[celltype$ClusterID %in% c(13),2]='Adipocyte' 

table(celltype$celltype)
sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.3
                      == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)

th <- theme(axis.text.x = element_text(angle = 60, hjust=1)) 
library(patchwork)
p_all_markers <- DotPlot(sce, features = genes_to_check, assay='RNA', 
                         group.by = 'celltype',  cols = c("#F7FBFF", "#08519C")) + coord_flip() + th
p_all_markers
ggsave(p_all_markers, filename = "./markers in celltype.pdf", width = 8, height = 10.5)

mycol <- colorRampPalette(brewer.pal(12, "Paired"))(11)
p_umap <- DimPlot(sce, reduction = "umap", group.by = "celltype", pt.size = 1.5,
                  cols = mycol, label = F, label.box = F)
p_umap
ggsave('./umap celltype 1.pdf', width = 6.5, height = 5)

library(plot1cell)
library(ggsci)
library(scales)
vignette("ggsci")
grep("^scale_color_", ls("package:ggsci"), value = TRUE)
show_col(pal_npg()(10))
# 1. JAMA
JAMA_11 <- colorRampPalette(ggsci::pal_jama()(7))(11)
show_col(JAMA_11)
# 2. AAAS
AAAS_11 <- colorRampPalette(ggsci::pal_aaas()(10))(11)
show_col(AAAS_11)
# 3. JCO
JCO_11 <- colorRampPalette(ggsci::pal_jco()(10))(11)
show_col(JCO_11)
# 4. NPG / Nature
NPG_11 <- colorRampPalette(ggsci::pal_npg()(10))(11)
show_col(NPG_11)
# 5. Lancet
LANCET_11 <- colorRampPalette(ggsci::pal_lancet()(9))(11)
show_col(LANCET_11)
# 6. NEJM
NEJM_11 <- colorRampPalette(ggsci::pal_nejm()(8))(11)
show_col(NEJM_11)

colnames(sce@meta.data)  
###Prepare data for ploting
circ_data <- prepare_circlize_data(sce, scale = 0.8)
set.seed(1234)
cluster_colors <- rand_color(length(levels(sce)))

cluster_colors <- c("#1B1919B2", "#008B45B2", "#BB0021B2", "#A20056B2", "#808180B2",
                    "#3B4992B2", "#5F559BB2", "#631879B2", "#EFC000B2", "#EE0000B2", "#008280B2")

# group_colors <- rand_color(length(names(table(sce$group))))
# group_colors <- c("#3C5488", "#E64B35", "#F28E2B", "#00A087")
rep_colors <- rand_color(length(names(table(sce$orig.ident))))
# rep_colors <- c("#AC1321", "#1A1674", "#A76446", "#6D6D6D")  # "#39004B",
rep_colors <- c("#FFD966B2", "#6FA8DCB2", "#8E7CC3B2", "#E06666B2")
orig_levels <- names(table(circ_data$orig.ident))
data.frame(orig.ident = orig_levels, color = rep_colors)

###plot and save figures
png(filename =  'circlize_plot.png', width = 6, height = 6,units = 'in', res = 300)
pdf("./umap_circle_1.pdf", width = 5, height = 5)
plot_circlize(circ_data, do.label = T, pt.size = 0.1, col.use = cluster_colors, 
              bg.color = "white", kde2d.n = 300, contour.nlevels = 100, repel = T, label.cex = 0.6)
# add_track(circ_data, group = "Group", colors = group_colors, tCELL_11# add_track(circ_data, group = "Group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "orig.ident", colors = rep_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
dev.off()

# 某个基因的表达情况
FeaturePlot(sce.all,features = c('OPA1', 'YY1', 'CEBPB' ,'NRBF2'), 
            ncol = 2, label.size = 3)

VlnPlot(test, features = c('OPA1', 'YY1', 'CEBPB', 'SLC5A2'),
        group.by = 'celltype', pt.size=0)
# geom_boxplot()+NoLegend()

################## step8:细胞数目和比例对比 ################## 
### 找各个单细胞亚群特异性高表达量基因
sce <- sce.all
table(Idents(sce))

sample.markers <- FindAllMarkers(sce, only.pos = FALSE, 
                                 min.pct = 0.25, logfc.threshold = 0.25)
top5 <- sample.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_table <- unstack(top5, gene ~ cluster)
jjVolcano(diffData = sample.markers,
          log2FC.cutoff = 0.5,
          adjustP.cutoff = 0.5,
          tile.col = celltype_colors,
          # tile.col = endo_sub_col,
          back.col='white',
          pSize = 0.4,
          aesCol = c("#0099CCBB","#CC3333BB"),
          # aesCol = c('purple','orange'),
          # aesCol = c("#8BB399", "#E1BB74"),
          # aesCol = c("#676C44",  "#E7CEC7"),
          celltypeSize = 3,
          # legend.position=c(0.1,0.9)
          topGeneN = 3)  #+ ylim(-5,8)
ggsave(filename = "./volcano_sub_celltype_top5gene.pdf", width = 11.5, height = 8)

library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(ggtext)

# 构建统计数据
cell_df <- data.frame(
  group = as.character(sce$sample),          
  celltype = as.character(Idents(sce)),
  stringsAsFactors = FALSE)

cell_df$group <- factor(cell_df$group, levels = c("50", "100", "150", "200"))

cell_stat <- cell_df %>%
  dplyr::count(group, celltype, name = "n") %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    proportion = n / sum(n),
    percent = proportion * 100
  ) %>%
  dplyr::ungroup()

# 细胞类型排序：这里按平均比例排序（从高到低）
cell_order <- cell_stat %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(mean_percent = mean(percent), .groups = "drop") %>%
  dplyr::arrange(desc(mean_percent)) %>%
  dplyr::pull(celltype)

# ggplot y轴从下往上画，因此这里反转
cell_levels <- rev(cell_order)
cell_stat$celltype <- factor(cell_stat$celltype, levels = cell_levels)

# group颜色（用于 50 / 100 / 150 / 200）
group_colors <- sample_colors

# 细胞类型颜色
cluster_colors2 <- cluster_colors

# 优化cell type名称显示
pretty_celltype <- c(
  "Adipocyte" = "Adipocyte",
  "B_cell" = "B cell",
  "Endothelial" = "Endothelial",
  "Fibroblast" = "Fibroblast",
  "Lymphatic_EC" = "Lymphatic EC",
  "Macrophage" = "Macrophage",
  "Monocyte" = "Monocyte",
  "Neutrophil" = "Neutrophil",
  "Proliferating_cell" = "Proliferating cell",
  "Smooth_muscle_cell" = "Smooth muscle cell",
  "T_cell" = "T cell"
)

# 如果你的cell type还有别的名字，这里自动兜底
pretty_labels <- pretty_celltype[cell_levels]
pretty_labels[is.na(pretty_labels)] <- cell_levels[is.na(pretty_labels)]

# 构建带彩色圆点的markdown y轴标签
cell_labels_md <- setNames(
  paste0(
    "<span style='color:", cluster_colors2[cell_levels], "; font-size:25pt;'>●</span> ",
    "<span style='color:#222222;'>", pretty_labels, "</span>"
  ),
  cell_levels
) 

panel_theme_base <- theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(fill = NA, color = "#4D4D4D", linewidth = 1.15),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "#333333"),
    axis.title.x = element_text(size = 12, face = "bold", color = "#333333"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, color = "#222222"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "#222222"),
    plot.margin = margin(8, 8, 8, 8)
  )

# 左图：Cell number（向左）
p_count_final <- ggplot(cell_stat, aes(x = -n, y = celltype, fill = group)) +
  geom_col(position = position_dodge(width = 0.72),
    width = 0.62, color = NA, alpha = 0.95) +
  geom_vline(xintercept = 0, color = "#BDBDBD", linewidth = 0.5) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(labels = function(x) comma(abs(x)), expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_discrete(labels = cell_labels_md) +
  labs(x = "Cell number", y = NULL, title = "Count") +
  panel_theme_base +
  theme(panel.background = element_rect(fill = "#F3FAFC", color = NA),
    axis.text.y = ggtext::element_markdown(size = 11, lineheight = 1.05),
    legend.position = "top")

# 右图：Cell proportion（向右）
p_prop_final <- ggplot(cell_stat, aes(x = percent, y = celltype, fill = group)) +
  geom_col(position = position_dodge(width = 0.72),
    width = 0.62, color = NA, alpha = 0.95) +
  geom_vline(xintercept = 0, color = "#BDBDBD", linewidth = 0.5) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.03, 0.05))) +
  scale_y_discrete(labels = NULL) +
  labs(x = "Cell proportion", y = NULL, title = "Proportion") +
  panel_theme_base +
  theme(panel.background = element_rect(fill = "#FDF8FA", color = NA),
    axis.text.y = element_blank(),
    legend.position = "none")

p_final_submission <- p_count_final + p_prop_final +
  plot_layout(widths = c(1.15, 1), guides = "collect") &
  theme(legend.position = "top")

p_final_submission
ggsave("Celltype_count_proportion_final.pdf", p_final_submission, device = cairo_pdf, width = 7, height = 5.2, bg = "white")

#### Pieplot
data <- read.csv("./Pie_group_cell_counts_freq.csv")
sum <- data[,c(1,2)]
colnames(sum) <- c('Cell_type', 'Sum')
#构建数据
library(ggplot2)
library(ggforce)
#基本空心饼状图
pdf("./sum of celltype.pdf", width = 6.5, height = 4.7)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')+#添加色
  scale_fill_manual(values = celltype_colors)+
  #r0=1则为空心饼图，r0=0则为不空心饼图
  geom_arc_bar(data = sum, stat = "pie",
               aes(x0 = 0, y0 = 0, r0 = 1, r = 2, amount = Sum, fill = Cell_type),
               color = "white",size=1)   #改变边界线的粗细与颜色
dev.off()

######################## 各个细胞类型特异性高表达基因 ###########################
sample.markers <- FindAllMarkers(sce, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# 保存每个cluster top10的marker基因
top10_table=unstack(top10, gene ~ cluster)
names(top10_table)=gsub("X","cluster", names(top10_table))
write.csv(file="top10_marker_genes.csv", top10_table, row.names=F)
