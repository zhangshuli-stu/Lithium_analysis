############################################################################
## 单核 -> 巨噬转换：人鼠对照、调控因子、人类按样本指标
##   Part 1  S11 基因组在人鼠轨迹上的动态（人类版 S11 c/d，与小鼠并排）
##   Part 2  四个转换程序（单核 / 过渡 / 脂质重塑 / 铁处理）沿轨迹的时序，及随 LiCl 剂量的变化
##   Part 3  转录因子表达沿轨迹（离线：GO:0003700 转录因子名单），找人鼠一致的转录因子
##   Part 4  Wnt/β-catenin（锂的靶点 GSK3β 下游）沿轨迹与随剂量的变化
##   Part 5  人类按样本：Classical Mono 占比、程序得分，AAA vs Control
############################################################################
mouse_rds <- "/home/data/t080553/Li_AAA/03_myeloid/mono_mac_5state.rds"
cs_dir    <- "/home/data/t080553/Li_AAA/04_cross_species/results"
out_dir   <- "/home/data/t080553/Li_AAA/05_transition"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mouse_state_col <- "myeloid_subtype_short"
mouse_pt_col    <- "monocle_pt"
mouse_dose_col  <- "sample"                    # 剂量 50/100/150/200
mouse_gmt       <- "/home/data/t080553/Li_AAA/05_transition/mh.all.v2023.1.Mm.symbols.gmt"   # 小鼠 hallmark（Wnt 用）
manual_orth     <- c(Ccl8 = "CCL8", Ftl1 = "FTL")  # 非一对一但对应明确的同源基因，手动补充
mouse_lineage   <- c("Classical Mono", "Transitional Mac", "Remodeling Mac")
human_lineage   <- c("Classical_Mono", "Transitional_Mono", "TREM2_Foamy_Mac")
human_pt_col    <- "pt_sling_Remodeling"

## S11 的基因组（小鼠名），按程序分组
programs_mouse <- list(
  Monocyte     = c("Ccr2", "Ly6c2", "Plac8", "Sell"),
  Transitional = c("Gdf15", "Ccl8", "Adgre1", "Mrc1"),
  Remodeling   = c("Spp1", "Mmp12", "Trem2", "Gpnmb", "Lgals3", "Cd36", "Fabp5", "Arg1", "Igf1"),
  Iron         = c("Ftl1", "Fth1", "Slc40a1", "Hmox1"))
prog_cols <- c(Monocyte = "#D0915B", Transitional = "#7FB8C9", Remodeling = "#5B5FA6", Iron = "#B5543A")
dose_cols <- c("50" = "#FFD966", "100" = "#6FA8DC", "150" = "#8E7CC3", "200" = "#E06666")
species_cols <- c(Mouse = "#E69F00", Human = "#0072B2")
grp_cols <- c(Control = "#6FA8DC", AAA = "#E06666")
human_cols <- c("Classical_Mono" = "#D0915B", "Transitional_Mono" = "#B7DCE5", "Resident_LYVE1_Mac" = "#8DCBEB",
                "TREM2_Foamy_Mac" = "#5B5FA6", "TIMD4_CXCL13_Mac" = "#141A4F", "IFN_Mac" = "#8A74B6",
                "Inflammatory_Mac" = "#C8615D")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(ComplexHeatmap); library(circlize)
})
count <- dplyr::count; filter <- dplyr::filter; select <- dplyr::select; mutate <- dplyr::mutate
summarise <- dplyr::summarise; group_by <- dplyr::group_by; ungroup <- dplyr::ungroup; n <- dplyr::n
set.seed(123)
pdf_save <- function(p, name, w, h) ggsave(file.path(out_dir, name), p, width = w, height = h, bg = "white", limitsize = FALSE)
get_counts <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) { obj <- JoinLayers(obj); LayerData(obj, layer = "counts") }
  else GetAssayData(obj, slot = "counts")
}
th <- theme_classic(base_size = 10) + theme(axis.text = element_text(color = "black"),
                                            strip.background = element_blank(), strip.text = element_text(face = "bold"))

######################################################################
## 读取：小鼠（换成人类同源基因名）+ 人类
######################################################################
orth  <- read.csv(file.path(cs_dir, "orthologs_1to1.csv"))
mouse <- readRDS(mouse_rds)
mouse <- tryCatch(UpdateSeuratObject(mouse), error = function(e) mouse)
mcnt  <- get_counts(mouse)
orth  <- orth[orth$mouse %in% rownames(mcnt), ]
mcnt_h <- mcnt[orth$mouse, ]; rownames(mcnt_h) <- orth$human
extra <- intersect(names(manual_orth), rownames(mcnt))
if (length(extra)) {
  add <- mcnt[extra, , drop = FALSE]; rownames(add) <- manual_orth[extra]
  add <- add[!rownames(add) %in% rownames(mcnt_h), , drop = FALSE]
  mcnt_h <- rbind(mcnt_h, add)
}
mh <- CreateSeuratObject(mcnt_h, meta.data = data.frame(
  state = as.character(mouse@meta.data[[mouse_state_col]]), pt = mouse@meta.data[[mouse_pt_col]],
  dose = as.character(mouse@meta.data[[mouse_dose_col]]), row.names = colnames(mouse)))
mh <- NormalizeData(mh, verbose = FALSE)
mh$dose <- factor(mh$dose, names(dose_cols))
rm(mouse, mcnt, mcnt_h); gc()

ht <- readRDS(file.path(cs_dir, "C_human_trajectory_object.rds"))
hu <- readRDS(file.path(cs_dir, "human_with_cross_species_annotations.rds"))

## 程序基因：小鼠名 -> 人类名
to_h <- function(g) {
  h <- orth$human[match(g, orth$mouse)]
  miss <- is.na(h) & g %in% names(manual_orth)
  h[miss] <- manual_orth[g[miss]]
  h[!is.na(h)]
}
programs <- lapply(programs_mouse, to_h)
dropped <- setdiff(unlist(programs_mouse), c(orth$mouse, names(manual_orth)))
if (length(dropped)) message("无一对一人类同源基因，已去除：", paste(dropped, collapse = ", "))
panel_genes <- unique(unlist(programs))
gene_prog <- setNames(rep(names(programs), lengths(programs)), unlist(programs))

######################################################################
## Part 1：S11 基因组在人鼠轨迹上的动态
######################################################################
message("Part 1 ...")
m_lin <- mh[, mh$state %in% mouse_lineage & is.finite(mh$pt)]
h_lin <- ht[, ht$celltype %in% human_lineage & is.finite(ht@meta.data[[human_pt_col]])]
h_lin$pt <- h_lin@meta.data[[human_pt_col]]

bin_mean <- function(obj, pt, genes, nbin = 40) {
  b <- cut(rank(pt, ties.method = "first"), nbin, labels = FALSE)
  dat <- GetAssayData(obj, layer = "data")[genes, , drop = FALSE]
  m <- sapply(seq_len(nbin), function(i) Matrix::rowMeans(dat[, b == i, drop = FALSE]))
  rownames(m) <- genes; m
}
smooth01 <- function(m) {
  m <- t(apply(m, 1, function(v) { s <- stats::filter(v, rep(1/3, 3), sides = 2); s[1] <- s[2]; s[length(s)] <- s[length(s) - 1]; s }))
  t(apply(m, 1, function(v) (v - min(v)) / (max(v) - min(v) + 1e-9)))
}
g1 <- intersect(panel_genes, intersect(rownames(m_lin), rownames(h_lin)))
pm <- smooth01(bin_mean(m_lin, m_lin$pt, g1)); ph <- smooth01(bin_mean(h_lin, h_lin$pt, g1))
rownames(pm) <- rownames(ph) <- g1
g1 <- g1[order(match(gene_prog[g1], names(programs)))]
pm <- pm[g1, ]; ph <- ph[g1, ]
r_panel <- sapply(g1, function(g) suppressWarnings(cor(pm[g, ], ph[g, ])))
write.csv(data.frame(gene = g1, program = gene_prog[g1], r_mouse_human = round(r_panel, 3)),
          file.path(out_dir, "P1_panel_gene_curve_correlation.csv"), row.names = FALSE)

hcol <- colorRamp2(c(0, 0.5, 1), c("#2B5C9E", "#F5F5C6", "#C0392B"))
ra <- rowAnnotation(Program = gene_prog[g1], col = list(Program = prog_cols), show_annotation_name = FALSE)
h1 <- Heatmap(pm, name = "Scaled\nexpression", col = hcol, cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = factor(gene_prog[g1], names(programs)), row_title = NULL, left_annotation = ra,
              row_names_side = "left", row_names_gp = gpar(fontface = "italic", fontsize = 8), column_title = "Mouse")
h2 <- Heatmap(ph, col = hcol, show_heatmap_legend = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
              row_names_gp = gpar(fontface = "italic", fontsize = 8),
              right_annotation = rowAnnotation(r = anno_barplot(pmax(r_panel, 0), width = unit(1.2, "cm"),
                                                                gp = gpar(fill = prog_cols[gene_prog[g1]], col = NA)),
                                               annotation_name_side = "top"),
              column_title = "Human")
pdf(file.path(out_dir, "P1_panel_heatmap_mouse_vs_human.pdf"), width = 7.5, height = 6)
draw(h1 + h2, column_title = "Pseudotime (relative)  \u2192", column_title_side = "bottom")
dev.off()

cur <- rbind(data.frame(gene = rep(g1, ncol(pm)), x = rep(seq_len(ncol(pm)) / ncol(pm), each = length(g1)), y = as.vector(pm), species = "Mouse"),
             data.frame(gene = rep(g1, ncol(ph)), x = rep(seq_len(ncol(ph)) / ncol(ph), each = length(g1)), y = as.vector(ph), species = "Human"))
cur$gene <- factor(cur$gene, g1)
p <- ggplot(cur, aes(x, y, color = species)) + geom_line(linewidth = 0.8) + facet_wrap(~ gene, ncol = 6) +
  scale_color_manual(values = species_cols, name = NULL) + th + theme(legend.position = "top") +
  labs(x = "Relative pseudotime", y = "Scaled expression")
pdf_save(p, "P1_panel_curves_mouse_vs_human.pdf", 12, 1.8 * ceiling(length(g1) / 6) + 1)

######################################################################
## Part 2：四个转换程序沿轨迹的时序 + 随剂量的变化
######################################################################
message("Part 2 ...")
score_prog <- function(obj) {
  pr <- lapply(programs, intersect, rownames(obj)); pr <- pr[lengths(pr) >= 2]
  obj <- AddModuleScore(obj, features = pr, name = "prog_", seed = 123)
  for (i in seq_along(pr)) obj[[names(pr)[i]]] <- obj[[paste0("prog_", i)]]
  obj
}
mh <- score_prog(mh); hu <- score_prog(hu); h_lin <- score_prog(h_lin)
m_lin <- mh[, mh$state %in% mouse_lineage & is.finite(mh$pt)]
pn <- intersect(names(programs), colnames(m_lin@meta.data))

prog_curve <- function(obj, species) {
  md <- obj@meta.data; md$x <- rank(md$pt) / nrow(md)
  bind_rows(lapply(pn, function(p) data.frame(species = species, program = p, x = md$x, y = md[[p]])))
}
pc <- rbind(prog_curve(m_lin, "Mouse"), prog_curve(h_lin, "Human"))
pc$program <- factor(pc$program, pn); pc$species <- factor(pc$species, c("Mouse", "Human"))
p <- ggplot(pc, aes(x, y, color = program)) + geom_smooth(method = "loess", se = FALSE, linewidth = 1.1, span = 0.4) +
  facet_wrap(~ species, nrow = 1) + scale_color_manual(values = prog_cols, name = NULL) + th +
  labs(x = "Relative pseudotime", y = "Program score")
pdf_save(p, "P2_program_timing_mouse_vs_human.pdf", 8, 3.5)

## 小鼠：拟时序分布随剂量
md <- m_lin@meta.data
p1 <- ggplot(md, aes(pt, color = dose, fill = dose)) + geom_density(alpha = 0.15, linewidth = 0.9) +
  scale_color_manual(values = dose_cols, name = "LiCl") + scale_fill_manual(values = dose_cols, name = "LiCl") +
  th + labs(x = "Pseudotime", y = "Density")
early <- md %>% group_by(dose) %>% summarise(early = mean(pt <= quantile(md$pt, 0.25)) * 100, .groups = "drop")
p2 <- ggplot(early, aes(dose, early, fill = dose)) + geom_col(width = 0.6) + scale_fill_manual(values = dose_cols, guide = "none") +
  th + labs(x = "LiCl (mg/kg)", y = "Cells in early pseudotime (%)")
pdf_save(p1 | p2, "P2_mouse_pseudotime_by_dose.pdf", 9, 3.5)
write.csv(early, file.path(out_dir, "P2_mouse_early_fraction_by_dose.csv"), row.names = FALSE)

## 小鼠：巨噬细胞中各程序得分随剂量
mac_m <- mh@meta.data[mh$state %in% c("Transitional Mac", "Remodeling Mac", "Resident Mac"), ]
pd <- bind_rows(lapply(pn, function(p) data.frame(program = p, dose = mac_m$dose, y = mac_m[[p]])))
pd$program <- factor(pd$program, pn)
p <- ggplot(pd, aes(dose, y, fill = dose)) + geom_violin(scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") + facet_wrap(~ program, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = dose_cols, guide = "none") + th + labs(x = "LiCl (mg/kg)", y = "Program score (macrophages)")
pdf_save(p, "P2_mouse_program_by_dose.pdf", 10, 3.2)

######################################################################
## Part 3：转录因子表达沿轨迹（离线版：GO:0003700 转录因子名单）
######################################################################
message("Part 3 ...")
tf_list <- unique(suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = "GO:0003700",
                                                         keytype = "GOALL", columns = "SYMBOL"))$SYMBOL)
tf_g  <- intersect(tf_list, intersect(rownames(m_lin), rownames(h_lin)))
det_m <- Matrix::rowMeans(GetAssayData(m_lin, layer = "data")[tf_g, ] > 0)
det_h <- Matrix::rowMeans(GetAssayData(h_lin, layer = "data")[tf_g, ] > 0)
tf_g  <- tf_g[det_m > 0.05 & det_h > 0.05]
cat("候选转录因子（两物种均有表达）：", length(tf_g), "\n")

tm  <- bin_mean(m_lin, m_lin$pt, tf_g, nbin = 30)
th_ <- bin_mean(h_lin, h_lin$pt, tf_g, nbin = 30)
sm  <- function(m) { x <- t(apply(m, 1, function(v) as.numeric(stats::filter(v, rep(1/3, 3), sides = 2))))
x[, 1] <- m[, 1]; x[, ncol(m)] <- m[, ncol(m)]; rownames(x) <- rownames(m); x }
tm_s <- sm(tm); th_s <- sm(th_)
tf_res <- data.frame(tf = tf_g,
                     r_mouse_human = sapply(tf_g, function(t) suppressWarnings(cor(tm_s[t, ], th_s[t, ]))),
                     trend_mouse   = sapply(tf_g, function(t) suppressWarnings(cor(seq_len(ncol(tm)), tm_s[t, ], method = "spearman"))),
                     trend_human   = sapply(tf_g, function(t) suppressWarnings(cor(seq_len(ncol(th_)), th_s[t, ], method = "spearman"))),
                     fc_mouse      = apply(tm_s, 1, function(v) diff(range(v))),
                     fc_human      = apply(th_s, 1, function(v) diff(range(v))))
tf_res[is.na(tf_res)] <- 0
tf_res$conserved <- tf_res$r_mouse_human > 0.5 & sign(tf_res$trend_mouse) == sign(tf_res$trend_human) &
  abs(tf_res$trend_mouse) > 0.5 & abs(tf_res$trend_human) > 0.5
tf_res$direction <- ifelse(tf_res$trend_mouse > 0, "Up along trajectory", "Down along trajectory")
tf_res <- tf_res[order(-tf_res$conserved, -tf_res$r_mouse_human), ]
write.csv(tf_res, file.path(out_dir, "P3_TF_expression_conservation.csv"), row.names = FALSE)
cat("人鼠一致的转录因子：", sum(tf_res$conserved), "\n")
print(head(tf_res[tf_res$conserved, c("tf", "r_mouse_human", "trend_mouse", "trend_human", "direction")], 30))

top <- tf_res$tf[tf_res$conserved]
top <- head(top[order(-pmin(tf_res$fc_mouse, tf_res$fc_human)[match(top, tf_res$tf)])], 30)
if (length(top) >= 3) {
  top <- top[order(apply(tm_s[top, , drop = FALSE], 1, which.max))]
  z <- function(m) t(scale(t(m)))
  col_tf <- colorRamp2(c(-2, 0, 2), c("#3B5B92", "white", "#B2182B"))
  hh1 <- Heatmap(z(tm_s[top, ]), name = "Expression\n(z)", col = col_tf, cluster_rows = FALSE, cluster_columns = FALSE,
                 row_names_side = "left", row_names_gp = gpar(fontface = "italic", fontsize = 8), column_title = "Mouse")
  hh2 <- Heatmap(z(th_s[top, ]), col = col_tf, show_heatmap_legend = FALSE, cluster_rows = FALSE, cluster_columns = FALSE,
                 row_names_gp = gpar(fontface = "italic", fontsize = 8), column_title = "Human")
  pdf(file.path(out_dir, "P3_conserved_TF_expression.pdf"), width = 6.5, height = 1.5 + 0.2 * length(top))
  draw(hh1 + hh2, column_title = "Pseudotime (relative)  \u2192", column_title_side = "bottom"); dev.off()
}

######################################################################
## Part 4：Wnt/β-catenin 沿轨迹与随剂量
######################################################################
message("Part 4 ...")
gmt <- clusterProfiler::read.gmt(mouse_gmt)
wnt_m <- unique(as.character(gmt$gene[grepl("WNT_BETA_CATENIN", gmt$term)]))
wnt <- to_h(wnt_m)
cat("Wnt 基因集：小鼠", length(wnt_m), "个 -> 人类同源", length(wnt), "个\n")
mh <- AddModuleScore(mh, features = list(intersect(wnt, rownames(mh))), name = "Wnt_", seed = 123)
h_lin <- AddModuleScore(h_lin, features = list(intersect(wnt, rownames(h_lin))), name = "Wnt_", seed = 123)
m_lin <- mh[, mh$state %in% mouse_lineage & is.finite(mh$pt)]
wd <- rbind(data.frame(species = "Mouse", x = rank(m_lin$pt) / ncol(m_lin), y = m_lin$Wnt_1, dose = as.character(m_lin$dose)),
            data.frame(species = "Human", x = rank(h_lin$pt) / ncol(h_lin), y = h_lin$Wnt_1, dose = NA))
wd$species <- factor(wd$species, c("Mouse", "Human"))
p1 <- ggplot(wd, aes(x, y, color = species)) + geom_smooth(method = "loess", se = TRUE, span = 0.4) +
  scale_color_manual(values = species_cols, name = NULL) + th + labs(x = "Relative pseudotime", y = "Wnt/\u03b2-catenin score")
p2 <- ggplot(subset(wd, species == "Mouse"), aes(x, y, color = dose)) + geom_smooth(method = "loess", se = FALSE, span = 0.4) +
  scale_color_manual(values = dose_cols, name = "LiCl") + th + labs(x = "Relative pseudotime", y = "Wnt/\u03b2-catenin score")
pdf_save(p1 | p2, "P4_Wnt_along_trajectory.pdf", 9, 3.5)

######################################################################
## Part 5：人类按样本的转换指标，AAA vs Control
######################################################################
message("Part 5 ...")
mono_mac <- c("Classical_Mono", "Transitional_Mono", "Resident_LYVE1_Mac", "TREM2_Foamy_Mac",
              "TIMD4_CXCL13_Mac", "IFN_Mac", "Inflammatory_Mac")
hm <- hu@meta.data[hu$celltype %in% mono_mac & hu$group %in% c("Control", "AAA"), ]
hm$is_mac <- !hm$celltype %in% c("Classical_Mono", "Transitional_Mono")
samp <- hm %>% group_by(dataset, sample, group) %>%
  summarise(n_cells = n(),
            ClassicalMono_pct = mean(celltype == "Classical_Mono") * 100,
            Resident_pct      = mean(celltype == "Resident_LYVE1_Mac") * 100,
            TREM2_pct         = mean(celltype == "TREM2_Foamy_Mac") * 100,
            Mono_to_Mac_ratio = sum(!is_mac) / max(sum(is_mac), 1),
            across(any_of(pn), ~ mean(.x[is_mac]), .names = "{.col}_score_mac"), .groups = "drop")
samp <- samp[samp$n_cells >= 50, ]
samp$group <- factor(samp$group, c("Control", "AAA"))
write.csv(samp, file.path(out_dir, "P5_human_per_sample_metrics.csv"), row.names = FALSE)
mets <- setdiff(colnames(samp), c("dataset", "sample", "group", "n_cells"))
long <- tidyr::pivot_longer(samp, all_of(mets), names_to = "metric", values_to = "value")
long$metric <- factor(long$metric, mets)
pv <- long %>% group_by(metric) %>%
  summarise(p = tryCatch(wilcox.test(value ~ group)$p.value, error = function(e) NA), .groups = "drop")
long <- left_join(long, pv, by = "metric")
long$lab <- sprintf("%s\nP = %.3f", long$metric, long$p)
p <- ggplot(long, aes(group, value, fill = group)) + geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_point(aes(shape = dataset), position = position_jitter(width = 0.12), size = 2) +
  facet_wrap(~ lab, scales = "free_y", nrow = 2) + scale_fill_manual(values = grp_cols, guide = "none") +
  th + labs(x = NULL, y = NULL, shape = "Dataset")
pdf_save(p, "P5_human_per_sample_AAA_vs_Control.pdf", 12, 6)
print(as.data.frame(pv))
message("完成，结果在：", out_dir)
