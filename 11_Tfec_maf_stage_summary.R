############################################################################
## A. 疾病分期汇总图：瘤颈 -> 扩张 AAA -> 大瘤体 -> 破裂（读取 bulk 分析结果）
## B. TFEC / MAF 在人鼠单核-巨噬细胞中的定位
##    B1 各亚群表达   B2 沿分化轨迹（人鼠对照）   B3 小鼠随 LiCl 剂量   B4 人类 AAA vs 对照（按样本）
## 全部离线运行
############################################################################
bulk_dir  <- "/home/data/t080553/Li_AAA/06_clinical_bulk"
cs_dir    <- "/home/data/t080553/Li_AAA/04_cross_species"
mouse_rds <- "/home/data/t080553/Li_AAA/03_myeloid/mono_mac_5state.rds"
out_dir   <- "/home/data/t080553/Li_AAA/07_TFEC_MAF"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork) })
filter <- dplyr::filter; mutate <- dplyr::mutate; summarise <- dplyr::summarise; group_by <- dplyr::group_by
select <- dplyr::select; ungroup <- dplyr::ungroup; n <- dplyr::n
th <- theme_classic(base_size = 10) + theme(axis.text = element_text(color = "black"),
                                            strip.background = element_blank(), strip.text = element_text(face = "bold"))
pdf_save <- function(p, name, w, h) ggsave(file.path(out_dir, name), p, width = w, height = h, bg = "white", limitsize = FALSE)

human_cols <- c("Classical_Mono" = "#D0915B", "Transitional_Mono" = "#B7DCE5", "Resident_LYVE1_Mac" = "#8DCBEB",
                "TREM2_Foamy_Mac" = "#5B5FA6", "TIMD4_CXCL13_Mac" = "#141A4F", "IFN_Mac" = "#8A74B6",
                "Inflammatory_Mac" = "#C8615D", "cDC2" = "#74A793", "cDC1" = "#3B7A5D")
mouse_cols <- c("Classical Mono" = "#D0915B", "Transitional Mac" = "#B7DCE5", "Resident Mac" = "#8DCBEB",
                "Remodeling Mac" = "#5B5FA6", "Osteoclast-like Mac" = "#141A4F")
dose_cols    <- c("50" = "#FFD966", "100" = "#6FA8DC", "150" = "#8E7CC3", "200" = "#E06666")
grp_cols     <- c(Control = "#6FA8DC", AAA = "#E06666")
species_cols <- c(Mouse = "#E69F00", Human = "#0072B2")
pretty <- function(x) gsub("_", " ", x)

######################################################################
## A. 疾病分期汇总图
######################################################################
message("A ...")
res  <- read.csv(file.path(bulk_dir, "comparisons_per_dataset.csv"))
meta <- read.csv(file.path(bulk_dir, "comparisons_meta.csv"))
stages <- c(Neck_vs_Control   = "Neck vs Control\n(GSE47472)",
            AAA_vs_Control    = "Dilated AAA vs Control\n(4 cohorts)",
            Large_vs_Small    = "Large vs Small AAA\n(2 cohorts)",
            Ruptured_vs_Large = "Ruptured vs Large AAA\n(GSE98278)")
ms <- c(Retention_index = "Monocyte retention index", Late_program = "Late lipid/iron programme",
        adj_TFEC = "TFEC", adj_MAF = "MAF", adj_Iron_handling = "Iron handling")
sa <- bind_rows(lapply(names(stages), function(tp) bind_rows(lapply(names(ms), function(m) {
  mm <- meta[meta$type == tp & meta$metric == m & is.finite(meta$k) & meta$k > 1, ]
  if (nrow(mm)) return(data.frame(type = tp, metric = m, g = mm$g, lo = mm$lo, hi = mm$hi, p = mm$p))
  r <- res[res$type == tp & res$metric == m, ]
  if (!nrow(r)) return(NULL)
  data.frame(type = tp, metric = m, g = r$g, lo = r$g - 1.96 * sqrt(r$v), hi = r$g + 1.96 * sqrt(r$v), p = r$p)
}))))
sa$stage  <- factor(stages[sa$type], rev(stages))
sa$metric <- factor(ms[sa$metric], ms)
sa$star   <- ifelse(sa$p < 0.001, "***", ifelse(sa$p < 0.01, "**", ifelse(sa$p < 0.05, "*", "")))
sa$dir    <- ifelse(sa$g > 0, "Up", "Down")
write.csv(sa, file.path(out_dir, "A_stage_summary.csv"), row.names = FALSE)

pA1 <- ggplot(sa, aes(g, stage, color = dir)) + geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y") + geom_point(size = 2.6) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") +
  scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC"), guide = "none") +
  th + labs(x = "Hedges' g (95% CI)", y = NULL)
pdf_save(pA1, "Main_A_stage_summary_forest.pdf", 15, 3.6)

pA2 <- ggplot(sa, aes(metric, stage, fill = g)) + geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f%s", g, ifelse(star == "", "", paste0("\n", star)))), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-2, 2),
                       oob = scales::squish, name = "Hedges' g") +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank(),
             axis.text.x = element_text(angle = 30, hjust = 1)) + labs(x = NULL, y = NULL)
pdf_save(pA2, "Main_A_stage_summary_heatmap.pdf", 7, 3.8)

######################################################################
## B. 单细胞
######################################################################
message("B：读取单细胞 ...")
h_genes <- c("TFEC", "MAF"); m_genes <- c("Tfec", "Maf")
hu <- readRDS(file.path(cs_dir, "myeloid_subtype.rds")); hu <- NormalizeData(hu, verbose = FALSE)
lv <- intersect(names(human_cols), unique(as.character(hu$celltype)))
hu$celltype <- factor(as.character(hu$celltype), lv)
ht <- readRDS(file.path(cs_dir, "results", "C_human_trajectory_object.rds"))

mouse <- readRDS(mouse_rds); mouse <- tryCatch(UpdateSeuratObject(mouse), error = function(e) mouse)
mouse <- NormalizeData(mouse, verbose = FALSE)
mouse$state <- as.character(mouse$myeloid_subtype_short)
mouse$dose  <- factor(as.character(mouse$sample), names(dose_cols))
mouse$pt    <- mouse$monocle_pt
mouse$state <- factor(mouse$state, intersect(names(mouse_cols), unique(mouse$state)))

## ---- B1 各亚群表达 ----
message("B1 ...")
d1 <- DotPlot(hu, features = h_genes, group.by = "celltype") + coord_flip() +
  scale_color_gradient(low = "grey90", high = "#B2182B") + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = human_cols[lv])) +
  scale_y_discrete(labels = pretty) + labs(x = NULL, y = NULL, title = "Human")
d2 <- DotPlot(mouse, features = m_genes, group.by = "state") + coord_flip() +
  scale_color_gradient(low = "grey90", high = "#B2182B") + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = mouse_cols[levels(mouse$state)])) +
  labs(x = NULL, y = NULL, title = "Mouse")
pdf_save(d1 | d2, "B1_TFEC_MAF_by_subtype_dotplot.pdf", 11, 3.6)
v1 <- VlnPlot(hu, features = h_genes, group.by = "celltype", cols = human_cols[lv], pt.size = 0, ncol = 2) &
  theme(axis.title.x = element_blank())
pdf_save(v1, "B1_human_TFEC_MAF_violin.pdf", 11, 4)

## ---- B2 沿分化轨迹（人鼠对照） ----
message("B2 ...")
h_lin <- ht[, ht$celltype %in% c("Classical_Mono", "Transitional_Mono", "TREM2_Foamy_Mac") &
              is.finite(ht$pt_sling_Remodeling)]
h_lin <- NormalizeData(h_lin, verbose = FALSE)
m_lin <- mouse[, mouse$state %in% c("Classical Mono", "Transitional Mac", "Remodeling Mac") & is.finite(mouse$pt)]
traj_df <- function(obj, genes, pt, species, label) {
  e <- FetchData(obj, vars = genes); colnames(e) <- label
  x <- rank(pt) / length(pt)
  bind_rows(lapply(label, function(g) data.frame(species = species, gene = g, x = x, y = e[[g]])))
}
td <- rbind(traj_df(h_lin, h_genes, h_lin$pt_sling_Remodeling, "Human", h_genes),
            traj_df(m_lin, m_genes, m_lin$pt, "Mouse", h_genes))
td$species <- factor(td$species, c("Mouse", "Human"))
pB2 <- ggplot(td, aes(x, y, color = species, fill = species)) + geom_smooth(method = "loess", span = 0.4, alpha = 0.2) +
  facet_wrap(~ gene, scales = "free_y") + scale_color_manual(values = species_cols, name = NULL) +
  scale_fill_manual(values = species_cols, name = NULL) + th +
  labs(x = "Relative pseudotime (Classical Mono \u2192 Remodeling/TREM2\u207a Mac)", y = "Expression")
pdf_save(pB2, "B2_TFEC_MAF_along_trajectory.pdf", 8, 3.5)

## ---- B3 小鼠随 LiCl 剂量 ----
message("B3 ...")
mac_states <- c("Transitional Mac", "Remodeling Mac", "Resident Mac")
mm <- FetchData(mouse, vars = c(m_genes, "dose", "state"))
mm <- mm[mm$state %in% mac_states, ]
dose_tab <- bind_rows(lapply(m_genes, function(g) mm %>% group_by(dose) %>%
                               summarise(gene = g, mean = mean(.data[[g]]), pct_pos = mean(.data[[g]] > 0) * 100, n = n(), .groups = "drop")))
write.csv(dose_tab, file.path(out_dir, "B3_mouse_TFEC_MAF_by_dose.csv"), row.names = FALSE)
cat("\n小鼠巨噬细胞中 Tfec / Maf 随剂量：\n"); print(as.data.frame(dose_tab))
long_m <- bind_rows(lapply(m_genes, function(g) data.frame(gene = g, dose = mm$dose, y = mm[[g]], state = mm$state)))
pB3a <- ggplot(long_m, aes(dose, y, fill = dose)) + geom_violin(scale = "width", linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  facet_wrap(~ gene, scales = "free_y") + scale_fill_manual(values = dose_cols, guide = "none") +
  th + labs(x = "LiCl (mg/kg)", y = "Expression (macrophages)")
pB3b <- ggplot(dose_tab, aes(dose, pct_pos, fill = dose)) + geom_col(width = 0.6) + facet_wrap(~ gene) +
  scale_fill_manual(values = dose_cols, guide = "none") + th + labs(x = "LiCl (mg/kg)", y = "% positive macrophages")
## 分亚群看：排除各剂量亚群组成不同的影响
by_state <- long_m %>% group_by(gene, state, dose) %>% summarise(mean = mean(y), n = n(), .groups = "drop") %>% filter(n >= 30)
pB3c <- ggplot(by_state, aes(dose, mean, color = state, group = state)) + geom_line(linewidth = 0.9) + geom_point(size = 2) +
  facet_wrap(~ gene, scales = "free_y") + scale_color_manual(values = mouse_cols, name = NULL) +
  th + labs(x = "LiCl (mg/kg)", y = "Mean expression")
pdf_save((pB3a | pB3b) / pB3c, "B3_mouse_TFEC_MAF_by_dose.pdf", 10, 7)
## 沿轨迹分剂量
td_m <- bind_rows(lapply(m_genes, function(g) {
  e <- FetchData(m_lin, vars = g)[[1]]
  data.frame(gene = g, dose = m_lin$dose, x = rank(m_lin$pt) / ncol(m_lin), y = e) }))
pB3d <- ggplot(td_m, aes(x, y, color = dose)) + geom_smooth(method = "loess", span = 0.4, se = FALSE) +
  facet_wrap(~ gene, scales = "free_y") + scale_color_manual(values = dose_cols, name = "LiCl") + th +
  labs(x = "Relative pseudotime", y = "Expression")
pdf_save(pB3d, "B3_mouse_TFEC_MAF_along_trajectory_by_dose.pdf", 8, 3.5)

## ---- B4 人类 AAA vs 对照（按样本，巨噬细胞） ----
message("B4 ...")
mac_h <- c("Resident_LYVE1_Mac", "TREM2_Foamy_Mac", "TIMD4_CXCL13_Mac", "IFN_Mac", "Inflammatory_Mac")
hh <- FetchData(hu, vars = c(h_genes, "celltype", "sample", "group", "dataset"))
hh <- hh[hh$celltype %in% mac_h & hh$group %in% c("Control", "AAA"), ]
samp <- hh %>% group_by(dataset, sample, group) %>%
  summarise(n_cells = n(), TFEC = mean(TFEC), MAF = mean(MAF), .groups = "drop") %>% filter(n_cells >= 30)
samp$group <- factor(samp$group, c("Control", "AAA"))
write.csv(samp, file.path(out_dir, "B4_human_TFEC_MAF_per_sample.csv"), row.names = FALSE)
lp <- bind_rows(lapply(h_genes, function(g) data.frame(gene = g, samp[, c("dataset", "sample", "group")], y = samp[[g]])))
pv <- lp %>% group_by(gene) %>% summarise(p = tryCatch(wilcox.test(y ~ group)$p.value, error = function(e) NA), .groups = "drop")
lp <- left_join(lp, pv, by = "gene"); lp$lab <- sprintf("%s (P = %.3f)", lp$gene, lp$p)
pB4 <- ggplot(lp, aes(group, y, fill = group)) + geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_point(aes(shape = dataset), position = position_jitter(width = 0.12), size = 2) +
  facet_wrap(~ lab, scales = "free_y") + scale_fill_manual(values = grp_cols, guide = "none") +
  th + labs(x = NULL, y = "Mean expression in macrophages (per sample)", shape = "Dataset")
pdf_save(pB4, "B4_human_TFEC_MAF_AAA_vs_Control.pdf", 7, 3.8)
print(as.data.frame(pv))
message("完成，结果在：", out_dir)
