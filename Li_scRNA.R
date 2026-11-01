.libPaths()
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

################## 配色选择 ################## 
# 配色
# 连续型sequential(连续的),：生成一系列连续渐变的颜色（共18组颜色，每组9个渐变色），
# 通常用来标记连续型数值的大小，适用于顺序数据。
display.brewer.all(type = "seq")
# 离散型diverging(离散的)：生成用深色强调两端、浅色标示中部的系列颜色（共9组颜色，每组11个颜色），
# 可用来标记数据中的离群点，适用于突显极端数值。
display.brewer.all(type = "qual")
# 极端型(Diverging)
display.brewer.all(type="div")
# 确定配色
mycol <- colorRampPalette(brewer.pal(12, "Paired"))(12)
cols = brewer.pal(9, name = "RdPu")
brewer.pal(9, name = "RdPu")
# "#FFF7F3" "#FDE0DD" "#FCC5C0" "#FA9FB5" "#F768A1" "#DD3497" "#AE017E" "#7A0177" "#49006A"
brewer.pal(9, name = "GnBu")
# "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"
# "#F7FCF0" "#E0F3DB" "#CCEBC5" "#A8DDB5" "#7BCCC4" "#4EB3D3" "#2B8CBE" "#0868AC" "#084081"

my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

my13colors <- c('#E5D2DD', '#D6E7A3', '#E0D4CA', '#91D0BE', '#CCE0F5',
                '#F1BB72', '#F3B1A0', '#6778AE', '#C5DEBA', '#CCC9E6',
                '#F7F398', '#C1E6F3', '#DCC1DD')

################## step1:导入数据 ##################  
getwd()
setwd("~/home/data/t080553/Li_AAA")

fs <- list.files('./scRNA_li_AAA/')
#执行这一步需要解压tar -xvf
samples <- str_split(fs,'in',simplify = T)[,1]


dir <- '~/Documents/R-training/Li/scRNA_li_AAA/'
# 每个样本创建一个Seurat对象
sceList = lapply(fs,function(pro){ 
  print(pro)  
  sce = CreateSeuratObject(counts =  Read10X(file.path(dir, pro)) ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)
})

# 31894 features across 16051
# 31918 features across 16033

sce.all <- merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = c("100_Li", "200_Li"))
# sce.all: 18245 x 12595
sce <- sce.all
saveRDS(sce, "./scRNA_li_AAA/00_RAW/sce_raw.rds")

sce <- readRDS("./00_RAW/sce_raw.rds")

as.data.frame(sce@assays$RNA@layers$counts.100[1:10, 1:2])
head(sce@meta.data, 10)
tail(sce@meta.data, 10)
table(sce@meta.data$orig.ident) 

sce$group <- ifelse(grepl("100",sce$orig.ident), "Li_100mg", "Li_200mg")
# sce$group[which(str_detect(sce$orig.ident, "GSM5077732"))] <- "Normal"
table(sce$group)

################## step2:质控 ##################  
setwd("01_QC/")
mito_genes <- rownames(sce)[grep("^mt-", rownames(sce))] 
mito_genes # 13线粒体基因
# 计算线粒体百分比
sce <- PercentageFeatureSet(sce, pattern = "^mt-", col.name = "percent_mito")
# sce.all[["percent.mito"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
fivenum(sce@meta.data$percent_mito)
# 检查数据的可用性
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), 
        group.by = "orig.ident",
        ncol = 3, cols = brewer.pal(9,"Set1")[1:6])
ggsave(filename = "./Vlnplot_feature_count_mito.pdf", width = 20, height = 8)

plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent_mito", 
                        group.by = "orig.ident", cols = brewer.pal(9,"Set1")[1:6])
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                        group.by = "orig.ident", cols = brewer.pal(9,"Set1")[1:6])
plot1 + plot2
ggsave(filename = "./FeatureScatter_feature_count_mito.pdf", width = 10, height = 5)

# >25% of reads from mitochondria genes
# removal of cells with <200 or >2,500 genes detected
# resulting in 19,807 cells
sce2 <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mito < 25
               & nCount_RNA <= 15000)
sce <- sce2
saveRDS(sce, "./sce_QC.rds")
sce <- readRDS("./sce_QC.rds")

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
ggsave(p1+p1_1, filename="./Dimplot_harmony_afterQC.pdf", width=10, height = 5)
dev.off()

################## step4: 降维聚类分群 ################## 
setwd("./01_QC/")
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
ggsave(p, filename = "./RNA_snn_res.0.2-0.8.pdf", width = 22, height = 10) # 分辨率为0.4,"RNA_snn_res.0.4"

# 查看特定分辨率的分群效果
colnames(sce@meta.data)
table(sce$RNA_snn_res.0.2)

# 删除一些RNA_snn_res
sce@meta.data <- sce@meta.data[ , -which(colnames(sce@meta.data) %in% c(
  "RNA_snn_res.0.3", "RNA_snn_res.0.4", "RNA_snn_res.0.5",
  "RNA_snn_res.0.6","RNA_snn_res.0.7", "RNA_snn_res.0.8"))]

# 最终降纬参数：dim=15, res=0.2
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:15)
sce <- FindClusters(sce, resolution = 0.2)
sce@reductions$umap <- NULL
sce@reductions$umap_naive <- NULL
sce@reductions$umap_harmony <- NULL
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:25, 
                 reduction.name = 'umap', reduction.key = 'UMAP_')
table(sce@meta.data$seurat_clusters)
sel.clust = "RNA_snn_res.0.2"
# sel.clust = "celltype"
sce <- SetIdent(sce, value = sel.clust)
table(sce@active.ident) 
ncluster <- length(unique(sce[[]]$RNA_snn_res.0.2))
mycol <- colorRampPalette(brewer.pal(12, "Paired"))(ncluster)
plot1 <- DimPlot(sce, reduction = "umap", cols = mycol, label = T, repel = F, 
                 pt.size = 2, shuffle = T)
ggsave(plot1, filename = "./umap.pdf", width = 4.5, height = 4) 
saveRDS(sce, "./sce.cluster.rds")

################## step5:简单检查常见分群情况 ################## 
setwd("../02_Celltype//")

DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.2", pt.size = 2, 
        cols = mycol, label = T) 

table(sce@meta.data$group)
DimPlot(sce, reduction = "umap", split.by = "group", cols = mycol, pt.size = 2, label = TRUE)
ggsave(filename = './UMAP_multi_samples_split.pdf', width = 9, height = 5)

# 查看细胞周期影响
FeaturePlot(sce, features = c("Pcna","Mki67"), reduction = "umap")
##s期基因PCNA，G2/M期基因MKI67均高表达，说明细胞均处于细胞周期
ggsave(file = "Cycling_Pcna_Mki67_expression.pdf", width=12, height = 6)

# 鉴定marker基因
obj.combined <- sce
sample.markers <- FindAllMarkers(obj.combined, only.pos = FALSE, min.pct = 0.05, logfc.threshold = 0.1)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# 保存每个cluster top10的marker基因
top10_table <- unstack(top10, gene ~ cluster)
names(top10_table) <- gsub("X","cluster", names(top10_table))
write.csv(file="top10_marker_genes.csv", top10_table, row.names=F)

genes_to_check <- c(# 'Cidec', 'Adipoq',  # Adipocyte
                   # 'Alas2', 'Bpgm', 'Snca', 'Hba-a1',    # Eryth, RBC
                   'Trbc2', 'Cd3e', 'Cd3d', 'Cd3g',        # T cell, 'Cd8a', 
                   'Cd79a', 'Cd79b', 'Ly6d','Ms4a1',   # B_cells
                   'Mbp', 'Gpm6b', 'Cnp', 'Prnp',     # Neuron  
                   'Cd52', 'Rac2', 'Coro1a', 'Ltb',   # Macro, T cell
                   'Cd68', 'C1qa', 'C1qb', 'Pf4', 'Lyz2', 'Mrc1',  'Cd52', 'Rac2', 'Coro1a', 
                   # Mono/Macro, 'Flt3', 'Stmn1',
                   'Cdh5', 'Fabp4', 'Pecam1', 'Cldn5', # EC, 'Cytl1',
                   'Acta2', 'Tagln', 'Myl9', 'Myh11',  # VSMC, 'Aqp1',
                   'Serpinf1', 'Dcn', 'Gsn', 'Lum', 'Col1a1','Col3a1'# Fibro,'Sfrp4','Pdgfra', 
)

# genes_to_check=str_to_upper(unique(genes_to_check)) #转换为大写
genes_to_check <- rev(unique(genes_to_check))
# th <- theme(axis.text.x = element_text(angle = 90, hjust=1))
p_all_markers <- DotPlot(sce, features = genes_to_check, assay='RNA',
                         cols = c("#F7FBFF", "#08519C")) + coord_flip() # + th 
p_all_markers
ggsave(p_all_markers, filename = "./markers to check.pdf", width = 8, height = 10)

################## step6:分群情况可视化 ################## 
obj.combined <- sce
p_umap <- DimPlot(obj.combined, reduction = "umap", group.by = "rec_SID6",
                  label = F, label.box = F)  # cols = my36colors, 
p_all_markers+p_umap
ggsave(p_umap, 'markers_umap_3.pdf',width = 15,height = 8)

DimPlot(sce.all, reduction = "umap",
        group.by = 'group',label = T) 
ggsave('group_umap.pdf',width = 7,height = 6)

################## step7:细胞亚群命名 ################## 

# 0,8: mono
# 1: macro
# 2: fibro
# 3: mano/macro
# 4: EC
# 5: T cell
# 6: B cell
# 7: VSMC
# 9: Neuron

celltype <- data.frame(ClusterID=0:9, celltype= 0:9) 
#定义细胞亚群 
celltype[celltype$ClusterID %in% c(0,8),2]='Monocyte'
celltype[celltype$ClusterID %in% c(1),2]='Macrophage'
celltype[celltype$ClusterID %in% c(2),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(3),2]="Mono/Macro"
celltype[celltype$ClusterID %in% c(4),2]='Endothelial'  
celltype[celltype$ClusterID %in% c(5),2]='T cell'
celltype[celltype$ClusterID %in% c(6),2]='B cell' 
celltype[celltype$ClusterID %in% c(7),2]='Smooth muscle cell' 
celltype[celltype$ClusterID %in% c(9),2]='Neuron cell'

table(celltype$celltype)
sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.2
                      == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)
saveRDS(sce, "./sce.celltype.rds")

th <- theme(axis.text.x = element_text(angle = 60, hjust=1)) 
library(patchwork)
p_all_markers <- DotPlot(sce, features = genes_to_check, assay='RNA', 
                         group.by = 'celltype',  cols = c("#F7FBFF", "#08519C")) + coord_flip() + th
p_all_markers
ggsave(p_all_markers, filename = "./markers in celltype.pdf", width = 7, height = 8.5)

mycol <- colorRampPalette(brewer.pal(12, "Paired"))(10)
p_umap <- DimPlot(sce, reduction = "umap", group.by = "celltype", pt.size = 1.5,
                  cols = mycol, label = F, label.box = F)
p_umap
ggsave('umap celltype.pdf', width = 6, height = 5)


### 整体细胞数
cell.num <- table(Idents(obj.combined))
cell.freq <- round(prop.table(table(Idents(obj.combined)))*100,2)
cell.combined <- rbind(cell.num, cell.freq)
write.csv(file="5_combined_cell_counts_freq.csv",cell.combined)
### 分组计算细胞数目和比例
cell.num.group <- table(Idents(obj.combined), obj.combined$group)
colnames(cell.num.group)<-paste0(colnames(cell.num.group),'_cell_counts')
cell.freq.group <- round(prop.table(table(Idents(obj.combined), obj.combined$group), margin = 2) *100,2)
colnames(cell.freq.group)<-paste0(colnames(cell.freq.group),'_cell_Freq')
cell.group <- cbind(cell.num.group, cell.freq.group)
write.csv(file="5_group_cell_counts_freq.csv",cell.group)
### 表格和UMAP展示细胞数目变化
pdf('5_UMAP_multi_samples_split_anno.pdf',width = 17,height = 5)
p<-DimPlot(obj.combined, reduction = "umap", split.by = "group",label=T,repel=T) + NoLegend()
tb <- tableGrob(cell.group)
plot_grid(p, tb,ncol=2,rel_widths=c(0.6,0.4))
dev.off() 

######################## 各个细胞类型特异性高表达基因 ###########################
sample.markers <- FindAllMarkers(sce, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- sample.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# 保存每个cluster top10的marker基因
top10_table=unstack(top10, gene ~ cluster)
names(top10_table)=gsub("X","cluster", names(top10_table))
write.csv(file="top10_marker_genes.csv", top10_table, row.names=F)
