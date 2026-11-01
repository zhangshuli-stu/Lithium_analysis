.libPaths()

rm(list=ls())
options(stringsAsFactors = F) 

getwd()
setwd("/home/data/t080553/ROE/li_myeloid/")

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
library(Seurat)
library(SeuratData)
library(SeuratDisk)

group_colors <- c(
  "50"  = "#FFD966B2",
  "100" = "#6FA8DCB2",
  "150" = "#8E7CC3B2",
  "200" = "#E06666B2"
)

major_colors <- c(
  "Macrophage"    = "#BC3C29", 
  "Monocyte"      = "#38618B",
  "Neutrophil"    = "#5A7A7C")

lineage_colors <- c(
  "Monocyte"       = "#D39450",
  "Macrophage"     = "#5B78B8",
  "Fibrocyte-like" = "#8C6BB1",
  "Neutrophil"     = "#C56860",
  "Dendritic cell" = "#4E927A",
  "QC"             = "#9A9A9A"
)

subtype_colors <- c(
  ## Monocyte
  "Classical Mono" = "#D98C4A",
  ## Macrophage
  "Transitional Mac" = "#AFD8E2", "Resident Mac" = "#7ED0F0", 
  "Remodeling Mac"   = "#5A61A8", "Osteoclast-like Mac" = "#111455",
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

###################### ROGUE: Purity ########################
# Calculate the ROGUE value of each putative cluster for each sample

library(ROGUE)
library(tibble)

scRNA <- readRDS("./myeloid.rds")
Idents(scRNA) <- scRNA$celltype_major
metadata <- scRNA@meta.data
expr <- as.matrix(LayerData(scRNA, assay = "RNA", layer = "counts"))
dim(expr)
expr <- matr.filter(expr, min.cells = 10, min.genes = 10) 
ent.res <- SE_fun(expr)
SEplot(ent.res)
ggsave("SEplot.pdf", SEplot(ent.res), width = 5, height = 4)
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value
rogue.res <- rogue(expr, labels = metadata$celltype,
                   samples = metadata$group, platform = "UMI")
rogue.boxplot(rogue.res) 
rogue.res.0.9 <- rogue(expr, labels = metadata$celltype, samples = metadata$orig.ident, 
                       platform = "UMI", span = 0.9)
rogue.boxplot(rogue.res.0.9)

cell_n_major <- table(metadata$group, metadata$celltype) %>%
  as.data.frame() %>%
  setNames(c("sample", "celltype", "n_cells")) %>%
  mutate(across(c(sample, celltype), as.character))

plot_df_major <- as.data.frame(rogue.res) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "celltype", values_to = "ROGUE") %>%
  filter(!is.na(ROGUE)) %>%
  left_join(cell_n_major, by = c("sample", "celltype")) %>%
  mutate(celltype = factor(celltype, levels = names(major_colors)))

p <- ggplot(plot_df_major, aes(x = celltype, y = ROGUE)) +
  geom_boxplot(aes(fill = celltype, color = celltype),
               alpha = 0.25, width = 0.6,
               outlier.shape = NA, linewidth = 0.8, fatten = 2.2) +
  geom_jitter(aes(fill = celltype, color = celltype, size = n_cells),
              width = 0.13, shape = 21, stroke = 1.0, alpha = 0.45) +
  scale_fill_manual(values = major_colors) +
  scale_color_manual(values = major_colors) +
  scale_size_continuous(range = c(1.5, 5.5), trans = "log10",
                        breaks = c(10, 100, 1000),
                        name = "Cells per sample") +
  scale_y_continuous(limits = c(0.5, 0.8),
                     breaks = seq(0.5, 0.8, 0.1),
                     expand = expansion(mult = c(0.04, 0.06))) +
  guides(fill = "none", color = "none",
         size = guide_legend(override.aes = list(color = "grey40",
                                                 fill = NA, stroke = 1.0))) +
  labs(x = NULL, y = "ROGUE") +
  theme_bw(base_size = 13, base_family = "sans") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid   = element_blank(),
    axis.text.x  = element_text(color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.ticks   = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.length = unit(0.18, "cm"),
    legend.position = c(0.45, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.key = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.margin = margin(4, 6, 4, 6),
    plot.margin = margin(15, 10, 25, 15))
p
ggsave("Ro_score.pdf", p, width = 4, height = 4)


rogue.res.sub <- rogue(expr, labels = metadata$myeloid_subtype_short,
                   samples = metadata$group, platform = "UMI")
rogue.boxplot(rogue.res.sub) 
cell_n <- table(metadata$group, metadata$myeloid_subtype_short) %>%
  as.data.frame() %>%
  setNames(c("sample", "subtype", "n_cells")) %>%
  mutate(across(c(sample, subtype), as.character))
plot_df <- as.data.frame(rogue.res.sub) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "subtype", values_to = "ROGUE") %>%
  filter(!is.na(ROGUE), subtype %in% names(sub_colors)) %>%
  left_join(cell_n, by = c("sample", "subtype")) %>%
  mutate(subtype = factor(subtype, levels = names(sub_colors)))

p <- ggplot(plot_df, aes(x = subtype, y = ROGUE)) +
  geom_boxplot(aes(fill = subtype, color = subtype),
               alpha = 0.25, width = 0.6,
               outlier.shape = NA, linewidth = 0.8, fatten = 2.2) +
  geom_jitter(aes(fill = subtype, color = subtype, size = n_cells),
              width = 0.13, shape = 21, stroke = 1.0, alpha = 0.45) +
  scale_fill_manual(values = sub_colors) +
  scale_color_manual(values = sub_colors) +
  scale_size_continuous(range = c(1.5, 5.5), trans = "log10",
                        breaks = c(10, 100, 1000),
                        name = "Cells per sample") +
  scale_y_continuous(limits = c(0.5, 0.92),
                     breaks = seq(0.5, 0.9, 0.1),
                     expand = expansion(mult = c(0.04, 0.06))) +
  guides(fill = "none", color = "none") +
  labs(x = NULL, y = "ROGUE") +
  theme_bw(base_size = 13, base_family = "sans") +
  theme_bw(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid   = element_blank(),
    axis.text.x  = element_text(color = "black", size = 12,
                                angle = 40, hjust = 1, vjust = 1),
    axis.text.y  = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.ticks   = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.length = unit(0.18, "cm"),
    legend.position = c(0.28, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.key = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.margin = margin(4, 6, 4, 6),
    plot.margin = margin(15, 10, 25, 15))
p
ggsave("ROGUE_10subtypes.pdf", p, width = 6, height = 4.5)


###################### Ro/e: 单细胞亚群的组织偏好 ########################

library(tidyverse)

# 列联表：行=亚群，列=组织
tab <- table(metadata$myeloid_subtype_short, metadata$group)
# 卡方检验拿期望值
chisq_res <- chisq.test(tab)
roe <- chisq_res$observed / chisq_res$expected
round(roe, 2)

library(tidyverse)

drop_sub <- c("Platelet-associated", "Hb-high")

roe_df <- as.data.frame(roe) %>%
  setNames(c("subtype", "tissue", "Roe")) %>%
  filter(!subtype %in% drop_sub) %>%
  mutate(
    tissue  = factor(tissue, levels = c("50", "100", "150", "200")),
    subtype = factor(subtype, levels = rev(setdiff(names(sub_colors), drop_sub))),
    label = case_when(
      Roe >= 1.5 ~ "+++",
      Roe >= 1.2 ~ "++",
      Roe >= 1.0 ~ "+",
      Roe >= 0.8 ~ "±",
      Roe >  0   ~ "−",
      TRUE       ~ ""))

p <- ggplot(roe_df, aes(x = tissue, y = subtype, fill = Roe)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = label, color = Roe > 1.4),
            size = 4, fontface = "bold") +
  scale_fill_gradientn(
    colors = c(
      "#7186B8", "#B2BEDD", "#FAF8FB", "#E7AFC0", "#B95F7E"),
    # colors = c("#7C9BB5", "#B8CBD9", "#F4F1EC", "#E8C4B8", "#C97B6B"),
    values = scales::rescale(c(0, 0.6, 1, 1.4, 2)),
    limits = c(0, 2), oob = scales::squish,
    breaks = c(0, 1, 2),
    name = "Ro/e") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey25"),
                     guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  coord_fixed(ratio = 0.55) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    panel.grid   = element_blank(),
    axis.text.x  = element_text(color = "black", size = 11, margin = margin(t = 4)),
    axis.text.y  = element_text(color = "black", size = 11, margin = margin(r = 4)),
    axis.ticks   = element_blank(),
    legend.position = "right",
    legend.key.width  = unit(0.3, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10))
p
ggsave("Roe_heatmap.pdf", p, width = 5, height = 5)
