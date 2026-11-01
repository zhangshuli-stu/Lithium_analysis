############################################################################
## bulk 临床队列验证：单核 -> 巨噬转换特征与 AAA、瘤体大小、破裂的关系
##   Step 1  读取 6 个 bulk 数据集，统一为 log2 表达 + 分组
##   Step 2  构建特征：保守状态特征、轨迹程序、铁处理、单核滞留指数、关键转录因子
##   Step 3  打分（数据集内 z-score 均值）；并做泛巨噬含量校正（残差）
##   Step 4  比较：AAA vs 对照（主分析 4 个扩张瘤体队列；含瘤颈 5 队列为敏感性分析）、瘤颈 vs 对照、大 vs 小、破裂
##   Step 5  作图：森林图、分组箱线图（Main_ 开头为主图，只含预先设定的主要指标）
##   主要指标：单核滞留指数、晚期程序、TFEC 和 MAF（巨噬含量校正）；FDR 只在主要指标间校正
##   Step 6  TFEC 敏感性分析：同时校正髓系和 B 细胞含量
##   Step 7  单细胞支撑数据：TFEC 细胞来源（人类全组织）、小鼠各巨噬亚群 Tfec 随剂量（结果缓存为 csv）
##   Step 8  临床部分主图（Figure_clinical_main.pdf，a–f）与补充图（Supp_ 开头）
## 全部离线运行
############################################################################
bulk_dir <- "/home/data/t080553/Li_AAA/rawdata"
cs_dir   <- "/home/data/t080553/Li_AAA/04_cross_species"
out_dir  <- "/home/data/t080553/Li_AAA/06_clinical_bulk"
human_sc_dir <- "/home/data/t080553/Li_AAA/human_sc"                         # 人类全组织单细胞（Step 7）
mouse_rds    <- "/home/data/t080553/Li_AAA/03_myeloid/mono_mac_5state.rds"  # 小鼠单核/巨噬（Step 7）
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork) })
filter <- dplyr::filter; mutate <- dplyr::mutate; summarise <- dplyr::summarise; group_by <- dplyr::group_by
count <- dplyr::count; ungroup <- dplyr::ungroup; n <- dplyr::n; select <- dplyr::select
th <- theme_classic(base_size = 10) + theme(axis.text = element_text(color = "black"),
                                            strip.background = element_blank(), strip.text = element_text(face = "bold"))
pdf_save <- function(p, name, w, h) ggsave(file.path(out_dir, name), p, width = w, height = h, bg = "white", limitsize = FALSE)

######################################################################
## Step 1：读取 bulk 数据
######################################################################
files <- list.files(bulk_dir, pattern = "^bulk_GSE\\d+\\.rds$", full.names = TRUE)
bulk <- lapply(files, function(f) {
  x  <- readRDS(f)
  ex <- as.matrix(x$exp); storage.mode(ex) <- "numeric"
  if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)          # 未 log 的数据（FPKM / MAS5）转 log2
  ex <- ex[rowSums(is.na(ex)) == 0, , drop = FALSE]
  ex <- ex[apply(ex, 1, var) > 0, , drop = FALSE]
  grp  <- tolower(as.character(x$group))
  size <- if (!is.null(x$size)) tolower(as.character(x$size)) else rep(NA, ncol(ex))
  list(id = sub("^bulk_(GSE\\d+)\\.rds$", "\\1", basename(f)), exp = ex, group = grp, size = size)
})
names(bulk) <- sapply(bulk, `[[`, "id")
for (b in bulk) {
  cat("\n==", b$id, "：", nrow(b$exp), "基因 x", ncol(b$exp), "样本；log2 范围",
      round(min(b$exp), 1), "~", round(max(b$exp), 1), "\n")
  print(table(group = b$group, size = b$size, useNA = "ifany"))
}

######################################################################
## Step 2：构建特征
######################################################################
## 2a 保守状态特征：小鼠亚群 marker 的人类同源基因 ∩ 在对应人类亚群中富集（亚群间 z > 0.5）
## 小鼠亚群 marker（A1 保存的结果）-> 人类同源基因，每个亚群取前 50 个
mk   <- read.csv(file.path(cs_dir, "results", "A1_mouse_state_markers.csv"))
orth <- read.csv(file.path(cs_dir, "results", "orthologs_1to1.csv"))
mouse_sigs <- lapply(split(mk, mk$cluster), function(x) {
  x <- x[x$p_val_adj < 0.05, ]; x <- x[order(-x$avg_log2FC), ]
  g <- orth$human[match(x$gene, orth$mouse)]; head(g[!is.na(g)], 50) })
hu <- readRDS(file.path(cs_dir, "myeloid_subtype.rds"))
hu <- NormalizeData(hu, verbose = FALSE)
hdat <- GetAssayData(hu, layer = "data")
ct <- as.character(hu$celltype)
avg <- sapply(sort(unique(ct)), function(s) Matrix::rowMeans(expm1(hdat[, ct == s, drop = FALSE])))
zg  <- t(scale(t(log1p(avg)))); zg[is.na(zg)] <- 0
pick <- function(mouse_state, human_sub, thr = 0.5) {
  g <- intersect(mouse_sigs[[mouse_state]], rownames(zg)); g[zg[g, human_sub] > thr]
}
sigs <- list(
  Classical_Mono = pick("Classical Mono", "Classical_Mono"),
  Resident_Mac   = pick("Resident Mac",   "Resident_LYVE1_Mac"),
  Remodeling_Mac = pick("Remodeling Mac", "TREM2_Foamy_Mac"))

## 2b 轨迹程序（D 部分：M1 早期、M3 晚期脂质/铁）+ 铁处理基因
dyn_f <- file.path(cs_dir, "results", c("D_dynamic_genes.csv", "D_dynamic_genes_peak.csv"))
dyn <- read.csv(dyn_f[file.exists(dyn_f)][1])
sigs$Early_program <- dyn$gene[dyn$module == "M1"]
sigs$Late_program  <- dyn$gene[dyn$module == "M3"]
sigs$Iron_handling <- c("FTL", "FTH1", "SLC40A1", "HMOX1", "CD163")

## 2c 泛巨噬含量（用于校正）
pan_myeloid <- c("CD68", "CSF1R", "C1QA", "C1QB", "C1QC", "AIF1", "TYROBP", "FCER1G", "CD14", "MS4A7")

## 2c' B 细胞含量（TFEC 敏感性分析用：组织中约 20% 的 TFEC 来自 B 细胞）
b_genes <- c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1")

## 2d 关键转录因子（人鼠一致）
key_tfs <- c("MAF", "TFEC", "BHLHE41", "NFE2L1", "EPAS1", "SMAD3")

## 2e 预先设定的主要指标（来自单细胞与跨物种分析；只在这些指标之间做 FDR 校正）
##    单核滞留指数是两个髓系程序之差，本身已抵消巨噬含量，用原始值；单个转录因子用巨噬含量校正后的值
## 2f 取材部位：GSE47472 为瘤颈（未扩张）组织，预先从 AAA vs 对照主分析中排除，单独分析并作为敏感性分析
neck_ids <- c("GSE47472")

primary <- c("Retention_index", "Late_program", "adj_TFEC", "adj_MAF")
primary_lab <- c(Retention_index = "Monocyte retention index", Late_program = "Late lipid/iron programme",
                 adj_TFEC = "TFEC (myeloid-adjusted)", adj_MAF = "MAF (myeloid-adjusted)")

cat("\n特征基因数：\n"); print(lengths(sigs))
saveRDS(list(sigs = sigs, pan_myeloid = pan_myeloid, key_tfs = key_tfs), file.path(out_dir, "signatures.rds"))
write.csv(data.frame(signature = rep(names(sigs), lengths(sigs)), gene = unlist(sigs)),
          file.path(out_dir, "signature_genes.csv"), row.names = FALSE)

######################################################################
## Step 3：打分
######################################################################
score_set <- function(ex, genes, min_genes = 3) {
  g <- intersect(genes, rownames(ex))
  if (length(g) < min_genes) return(rep(NA, ncol(ex)))
  z <- t(scale(t(ex[g, , drop = FALSE]))); colMeans(z, na.rm = TRUE)
}
scores <- lapply(bulk, function(b) {
  s <- sapply(sigs, function(g) score_set(b$exp, g))
  s <- cbind(s, Pan_myeloid = score_set(b$exp, pan_myeloid), Bcell = score_set(b$exp, b_genes))
  s <- cbind(s, Retention_index = s[, "Classical_Mono"] - s[, "Late_program"])
  for (tf in key_tfs) s <- cbind(s, setNames(data.frame(if (tf %in% rownames(b$exp)) as.numeric(scale(b$exp[tf, ])) else NA), tf))
  s <- as.data.frame(s)
  ## 校正：对泛巨噬得分回归后的残差（反映巨噬细胞内部的状态构成）
  for (v in setdiff(colnames(s), "Pan_myeloid")) {
    ok <- is.finite(s[[v]]) & is.finite(s$Pan_myeloid)
    s[[paste0("adj_", v)]] <- NA
    if (sum(ok) > 4) s[[paste0("adj_", v)]][ok] <- resid(lm(s[[v]][ok] ~ s$Pan_myeloid[ok]))
  }
  ## 同时校正髓系和 B 细胞含量的 TFEC（敏感性分析）
  s$adj2_TFEC <- NA
  ok <- is.finite(s$TFEC) & is.finite(s$Pan_myeloid) & is.finite(s$Bcell)
  if (sum(ok) > 5) s$adj2_TFEC[ok] <- resid(lm(s$TFEC[ok] ~ s$Pan_myeloid[ok] + s$Bcell[ok]))
  s$adj_Bcell <- NULL
  cbind(dataset = b$id, sample = colnames(b$exp), group = b$group, size = b$size, s)
})
sc_all <- bind_rows(scores)
write.csv(sc_all, file.path(out_dir, "bulk_scores_all_samples.csv"), row.names = FALSE)
metrics <- setdiff(colnames(sc_all), c("dataset", "sample", "group", "size"))

######################################################################
## Step 4：比较
######################################################################
hedges_g <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]; n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(c(g = NA, v = NA, p = NA, n1 = n1, n2 = n2))
  sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  g <- (1 - 3 / (4 * (n1 + n2) - 9)) * (mean(x) - mean(y)) / sp
  c(g = g, v = (n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2)),
    p = suppressWarnings(wilcox.test(x, y)$p.value), n1 = n1, n2 = n2)
}
meta_dl <- function(g, v) {                       # DerSimonian-Laird 随机效应
  ok <- is.finite(g) & is.finite(v); g <- g[ok]; v <- v[ok]; k <- length(g)
  if (k < 2) return(c(g = NA, lo = NA, hi = NA, p = NA, I2 = NA, k = k))
  w <- 1 / v; gf <- sum(w * g) / sum(w); Q <- sum(w * (g - gf)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (v + tau2); gr <- sum(ws * g) / sum(ws); se <- sqrt(1 / sum(ws))
  c(g = gr, lo = gr - 1.96 * se, hi = gr + 1.96 * se, p = 2 * pnorm(-abs(gr / se)),
    I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) else 0, k = k)
}

## 定义比较
is_ctrl <- function(x) grepl("control|normal|ctrl|donor", x)
cmps <- list()
for (id in names(bulk)) {
  d <- sc_all[sc_all$dataset == id, ]
  if (any(is_ctrl(d$group)) && any(!is_ctrl(d$group))) {
    main_type <- if (id %in% neck_ids) "Neck_vs_Control" else "AAA_vs_Control"   # 主分析只含扩张瘤体
    cmps[[length(cmps) + 1]] <- list(type = main_type, dataset = id, case = !is_ctrl(d$group), ctrl = is_ctrl(d$group))
    cmps[[length(cmps) + 1]] <- list(type = "AAA_vs_Control_all", dataset = id, case = !is_ctrl(d$group), ctrl = is_ctrl(d$group))
  }
}
d <- sc_all[sc_all$dataset == "GSE57691", ]
if (nrow(d)) cmps[[length(cmps) + 1]] <- list(type = "Large_vs_Small", dataset = "GSE57691",
                                              case = grepl("large", d$group), ctrl = grepl("small", d$group))
d <- sc_all[sc_all$dataset == "GSE98278", ]
if (nrow(d)) {
  lab98 <- paste(d$group, d$size)
  rup <- grepl("rupt", lab98); lrg <- grepl("large", lab98) & !rup; mid <- grepl("intermediate", lab98) & !rup
  cmps[[length(cmps) + 1]] <- list(type = "Large_vs_Small", dataset = "GSE98278", case = lrg, ctrl = mid)
  cmps[[length(cmps) + 1]] <- list(type = "Ruptured_vs_Intact", dataset = "GSE98278", case = rup, ctrl = !rup)
  cmps[[length(cmps) + 1]] <- list(type = "Ruptured_vs_Large", dataset = "GSE98278", case = rup, ctrl = lrg)
}
cat("\n比较列表：\n")
for (cm in cmps) cat(sprintf("  %-20s %-10s case n = %d, control n = %d\n", cm$type, cm$dataset, sum(cm$case), sum(cm$ctrl)))

res <- bind_rows(lapply(cmps, function(cm) {
  d <- sc_all[sc_all$dataset == cm$dataset, ]
  bind_rows(lapply(metrics, function(m) {
    h <- hedges_g(d[[m]][cm$case], d[[m]][cm$ctrl])
    data.frame(type = cm$type, dataset = cm$dataset, metric = m, g = h["g"], v = h["v"], p = h["p"],
               n_case = h["n1"], n_ctrl = h["n2"])
  }))
}))
meta <- res %>% group_by(type, metric) %>%
  summarise(as.data.frame(t(meta_dl(g, v))), .groups = "drop") %>% as.data.frame()
meta$fdr <- ave(meta$p, meta$type, FUN = function(p) p.adjust(p, "BH"))          # 全部指标（探索性）
meta$fdr_primary <- NA
for (tp in unique(meta$type)) {
  i <- meta$type == tp & meta$metric %in% primary
  meta$fdr_primary[i] <- p.adjust(meta$p[i], "BH")
}
res$fdr_primary <- NA                                                              # 单队列比较（破裂）
for (tp in unique(res$type)) {
  i <- res$type == tp & res$metric %in% primary
  res$fdr_primary[i] <- p.adjust(res$p[i], "BH")
}
write.csv(res,  file.path(out_dir, "comparisons_per_dataset.csv"), row.names = FALSE)
write.csv(meta, file.path(out_dir, "comparisons_meta.csv"), row.names = FALSE)

## 打印主要结果
main_m <- c("Classical_Mono", "Resident_Mac", "Remodeling_Mac", "Early_program", "Late_program",
            "Iron_handling", "Retention_index", key_tfs)
show <- function(tp) {
  cat("\n================", tp, "================\n")
  r <- res[res$type == tp & res$metric %in% c(main_m, paste0("adj_", main_m), "Pan_myeloid"), ]
  if (length(unique(r$dataset)) > 1) {
    mm <- meta[meta$type == tp & meta$metric %in% r$metric, c("metric", "g", "lo", "hi", "p", "fdr", "I2", "k")]
    mm[, 2:7] <- round(mm[, 2:7], 3); print(mm[order(mm$p), ], row.names = FALSE)
  } else {
    rr <- r[, c("dataset", "metric", "g", "p", "n_case", "n_ctrl")]; rr$g <- round(rr$g, 2); rr$p <- signif(rr$p, 2)
    print(rr[order(rr$p), ], row.names = FALSE)
  }
}
for (tp in unique(res$type)) show(tp)

## 主要指标汇总
cat("\n================ 主要指标（预先设定）================\n")
pm_tab <- bind_rows(
  meta[meta$metric %in% primary & meta$k > 1, c("type", "metric", "g", "lo", "hi", "p", "fdr_primary", "I2", "k")],
  res[res$metric %in% primary & res$type %in% c("Neck_vs_Control", "Ruptured_vs_Intact", "Ruptured_vs_Large"),
      c("type", "metric", "g", "p", "fdr_primary")] %>% mutate(k = 1))
pm_tab[, c("g", "lo", "hi", "I2")] <- round(pm_tab[, c("g", "lo", "hi", "I2")], 3)
pm_tab$p <- signif(pm_tab$p, 3); pm_tab$fdr_primary <- signif(pm_tab$fdr_primary, 3)
print(pm_tab, row.names = FALSE)
write.csv(pm_tab, file.path(out_dir, "primary_metrics_summary.csv"), row.names = FALSE)

## 敏感性分析：AAA vs 对照，逐一剔除队列
loo <- bind_rows(lapply(unique(res$dataset[res$type == "AAA_vs_Control"]), function(drop) {
  r <- res[res$type == "AAA_vs_Control" & res$dataset != drop & res$metric %in% primary, ]
  r %>% group_by(metric) %>% summarise(as.data.frame(t(meta_dl(g, v))), .groups = "drop") %>% mutate(dropped = drop)
}))
loo <- as.data.frame(loo)
write.csv(loo, file.path(out_dir, "sensitivity_leave_one_out_AAA_vs_Control.csv"), row.names = FALSE)
cat("\n== 逐一剔除队列后的合并结果（AAA vs 对照，4 个扩张瘤体队列）==\n")
print(loo %>% mutate(g = round(g, 2), p = signif(p, 2)) %>% select(dropped, metric, g, p, I2), row.names = FALSE)

######################################################################
## Step 5：作图
######################################################################
## 5a 森林图：每个特征在各数据集及合并后的效应量（Hedges g，正值 = 在 case 组更高）
case_lab <- c(AAA_vs_Control = "Higher in AAA \u2192", AAA_vs_Control_all = "Higher in AAA \u2192",
              Large_vs_Small = "Higher in large AAA \u2192")
forest <- function(tp, ms, file, w = 9, labs_map = NULL, nrow = 2, h = 5) {
  r <- res[res$type == tp & res$metric %in% ms, ]
  r$lo <- r$g - 1.96 * sqrt(r$v); r$hi <- r$g + 1.96 * sqrt(r$v); r$what <- r$dataset
  mm <- meta[meta$type == tp & meta$metric %in% ms, ]
  if (nrow(mm) && length(unique(r$dataset)) > 1)
    r <- bind_rows(r, data.frame(type = tp, metric = mm$metric, g = mm$g, lo = mm$lo, hi = mm$hi, what = "Pooled", p = mm$p))
  r$what <- factor(r$what, rev(c(sort(unique(setdiff(r$what, "Pooled"))), "Pooled")))
  r$metric <- factor(r$metric, ms)
  if (!is.null(labs_map)) levels(r$metric) <- ifelse(levels(r$metric) %in% names(labs_map), labs_map[levels(r$metric)], levels(r$metric))
  r$lab <- ifelse(r$what == "Pooled", sprintf("P = %.3g", r$p), "")
  ggplot(r, aes(g, what, color = what == "Pooled")) + geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y") + geom_point(aes(shape = what == "Pooled"), size = 2.2) +
    geom_text(aes(label = lab), hjust = -0.15, vjust = -0.8, size = 2.6, color = "black") +
    facet_wrap(~ metric, nrow = nrow) + scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "#B2182B"), guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") + th +
    labs(x = paste0("Hedges' g   (", case_lab[tp], ")"), y = NULL, title = gsub("_", " ", tp)) -> p
  pdf_save(p, file, w, h)
}
st <- c("Classical_Mono", "Resident_Mac", "Remodeling_Mac", "Early_program", "Late_program", "Iron_handling", "Retention_index")
## 主图：主要指标
forest("AAA_vs_Control", primary, "Main_F1_forest_AAA_vs_Control_primary.pdf", 12, primary_lab, nrow = 1, h = 3.2)
forest("Large_vs_Small", primary, "Main_F2_forest_Large_vs_Small_primary.pdf", 12, primary_lab, nrow = 1, h = 2.6)
## 补充图：全部指标
forest("AAA_vs_Control_all", primary, "Supp_F1_forest_AAA_vs_Control_incl_neck_primary.pdf", 12, primary_lab, nrow = 1, h = 3.4)
forest("AAA_vs_Control", st, "F1_forest_AAA_vs_Control_raw.pdf", 13)
forest("AAA_vs_Control", paste0("adj_", st), "F1_forest_AAA_vs_Control_adjusted.pdf", 13)
forest("AAA_vs_Control", key_tfs, "F1_forest_AAA_vs_Control_TFs.pdf", 11)
forest("Large_vs_Small", c(st, paste0("adj_", st)), "F2_forest_Large_vs_Small.pdf", 16)

## 5b GSE98278：中等 -> 大 -> 破裂
d <- sc_all[sc_all$dataset == "GSE98278", ]
if (nrow(d)) {
  lab98 <- paste(d$group, d$size)
  d$stage <- factor(ifelse(grepl("rupt", lab98), "Ruptured", ifelse(grepl("large", lab98), "Large", "Intermediate")),
                    c("Intermediate", "Large", "Ruptured"))
  ms <- c(st, paste0("adj_", st), key_tfs)
  long <- bind_rows(lapply(ms, function(m) data.frame(metric = m, stage = d$stage, value = d[[m]])))
  pv <- long %>% group_by(metric) %>% summarise(
    p_trend = suppressWarnings(cor.test(as.numeric(stage), value, method = "spearman")$p.value),
    p_rup_vs_large = suppressWarnings(wilcox.test(value[stage == "Ruptured"], value[stage == "Large"])$p.value), .groups = "drop")
  long <- left_join(long, pv, by = "metric")
  long$lab <- sprintf("%s\ntrend P=%.3f; R vs L P=%.3f", long$metric, long$p_trend, long$p_rup_vs_large)
  long$lab <- factor(long$lab, unique(long$lab[order(match(long$metric, ms))]))
  p <- ggplot(long, aes(stage, value, fill = stage)) + geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.12, size = 1.2) + facet_wrap(~ lab, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c(Intermediate = "#9ECAE1", Large = "#4292C6", Ruptured = "#B2182B"), guide = "none") +
    th + theme(strip.text = element_text(size = 7)) + labs(x = NULL, y = "Score")
  pdf_save(p, "F3_GSE98278_size_rupture.pdf", 14, 2.6 * ceiling(length(ms) / 5))
  ## 主图：主要指标 + 铁处理
  ms2 <- c(primary, "adj_Iron_handling")
  lab2 <- c(primary_lab, adj_Iron_handling = "Iron handling (myeloid-adjusted)")
  l2 <- long[long$metric %in% ms2, ]
  l2$lab2 <- factor(sprintf("%s\ntrend P = %.3f\nRuptured vs Large P = %.3f", lab2[l2$metric], l2$p_trend, l2$p_rup_vs_large))
  l2$lab2 <- factor(l2$lab2, unique(l2$lab2[order(match(l2$metric, ms2))]))
  p2 <- ggplot(l2, aes(stage, value, fill = stage)) + geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.12, size = 1.2) + facet_wrap(~ lab2, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(Intermediate = "#9ECAE1", Large = "#4292C6", Ruptured = "#B2182B"), guide = "none") +
    th + theme(strip.text = element_text(size = 8)) + labs(x = NULL, y = "Score")
  pdf_save(p2, "Main_F3_GSE98278_size_rupture_primary.pdf", 13, 3.6)
  write.csv(as.data.frame(pv), file.path(out_dir, "GSE98278_stage_tests.csv"), row.names = FALSE)
}

######################################################################
## Step 6：TFEC 敏感性分析（同时校正髓系和 B 细胞含量）
######################################################################
message("Step 6 ...")
tf_vers <- c(TFEC = "Unadjusted", adj_TFEC = "Myeloid-adjusted", adj2_TFEC = "Myeloid + B cell-adjusted")
sens <- res[res$metric %in% names(tf_vers), ]
sens$lo <- sens$g - 1.96 * sqrt(sens$v); sens$hi <- sens$g + 1.96 * sqrt(sens$v)
mm <- meta[meta$metric %in% names(tf_vers) & is.finite(meta$k) & meta$k > 1, ]
sens_sum <- bind_rows(
  sens[sens$type %in% c("Neck_vs_Control", "Ruptured_vs_Intact", "Ruptured_vs_Large"), c("type", "metric", "g", "lo", "hi", "p")],
  mm[, c("type", "metric", "g", "lo", "hi", "p")])
sens_sum$version <- factor(tf_vers[sens_sum$metric], tf_vers)
write.csv(sens_sum, file.path(out_dir, "sensitivity_TFEC_Bcell_adjustment.csv"), row.names = FALSE)
cat("\n== TFEC 敏感性分析 ==\n")
print(sens_sum %>% mutate(g = round(g, 2), p = signif(p, 2)) %>% select(type, version, g, p), row.names = FALSE)
tf_cor <- sc_all %>% group_by(dataset) %>%
  summarise(r_TFEC_myeloid = cor(TFEC, Pan_myeloid, use = "complete.obs"),
            r_TFEC_Bcell   = cor(TFEC, Bcell, use = "complete.obs"), .groups = "drop")
write.csv(tf_cor, file.path(out_dir, "TFEC_correlation_with_myeloid_Bcell.csv"), row.names = FALSE)
print(as.data.frame(tf_cor))

######################################################################
## Step 7：单细胞支撑数据（结果缓存为 csv，重跑时直接读取）
######################################################################
message("Step 7 ...")
src_csv  <- file.path(out_dir, "sc_TFEC_cell_source.csv")
dose_csv <- file.path(out_dir, "sc_mouse_Tfec_by_dose_state.csv")
lineage_map <- c(Monocyte = "Myeloid", Macrophage = "Myeloid", DC = "Myeloid", pDC = "Myeloid", Neutrophil = "Myeloid",
                 Mast_cell = "Myeloid", B_cell = "B / plasma", Plasma_cell = "B / plasma", T_cell = "T / NK", NK_cell = "T / NK",
                 Smooth_muscle_cell = "Stromal", Fibroblast = "Stromal", Endothelial = "Stromal", Lymphatic_EC = "Stromal",
                 Pericyte = "Stromal", Adipocyte = "Stromal", Proliferating_cell = "Other")
if (!file.exists(src_csv)) {
  sc_files <- list.files(human_sc_dir, pattern = "_sce_celltype\\.rds$", full.names = TRUE)
  if (length(sc_files)) {
    dp <- list(); src <- list()
    ord_ct <- names(lineage_map)
    for (f in sc_files) {
      gse <- sub("_sce_celltype\\.rds$", "", basename(f)); obj <- readRDS(f)
      e <- FetchData(obj, vars = c("TFEC", "celltype"), layer = "data")
      src[[gse]] <- e %>% group_by(celltype) %>% summarise(TFEC = sum(expm1(TFEC)), .groups = "drop") %>%
        mutate(pct = TFEC / sum(TFEC) * 100, dataset = gse)
      obj$celltype <- factor(as.character(obj$celltype), rev(intersect(ord_ct, unique(as.character(obj$celltype)))))
      dp[[gse]] <- DotPlot(obj, features = c("TFEC", "CD68", "CD14", "MS4A1", "CD3E", "MYH11", "COL1A1", "PECAM1"),
                           group.by = "celltype") + scale_color_gradient(low = "grey90", high = "#B2182B") +
        scale_y_discrete(labels = function(x) gsub("_", " ", x)) + ggtitle(gse) + th +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) + labs(x = NULL, y = NULL)
      rm(obj); gc()
    }
    write.csv(bind_rows(src), src_csv, row.names = FALSE)
    pdf_save(wrap_plots(dp, nrow = 1), "Supp_TFEC_whole_tissue_dotplot.pdf", 16, 5)
  } else message("没找到人类全组织单细胞文件（", human_sc_dir, "），跳过 TFEC 细胞来源")
}
if (!file.exists(dose_csv) && file.exists(mouse_rds)) {
  mo <- readRDS(mouse_rds); mo <- tryCatch(UpdateSeuratObject(mo), error = function(e) mo)
  mo <- NormalizeData(mo, verbose = FALSE)
  md <- data.frame(Tfec = FetchData(mo, vars = "Tfec")[[1]], state = as.character(mo$myeloid_subtype_short),
                   dose = as.character(mo$sample))
  rm(mo); gc()
  dose_tab <- md %>% filter(state %in% c("Transitional Mac", "Remodeling Mac", "Resident Mac")) %>%
    group_by(state, dose) %>% summarise(mean = mean(Tfec), se = sd(Tfec) / sqrt(n()), pct_pos = mean(Tfec > 0) * 100,
                                        n = n(), .groups = "drop")
  write.csv(dose_tab, dose_csv, row.names = FALSE)
}

######################################################################
## Step 8：临床部分主图与补充图（统一版式）
######################################################################
message("Step 8 ...")
col_up <- "#B2182B"; col_dn <- "#2166AC"
stage_cols <- c(Intermediate = "#9ECAE1", Large = "#4292C6", Ruptured = "#B2182B")
dose_cols  <- c("50" = "#FFD966", "100" = "#6FA8DC", "150" = "#8E7CC3", "200" = "#E06666")
mstate_cols <- c("Transitional Mac" = "#7FB8C9", "Resident Mac" = "#8DCBEB", "Remodeling Mac" = "#5B5FA6")
lin_cols <- c(Myeloid = "#C8615D", "B / plasma" = "#74A793", "T / NK" = "#8A74B6", Stromal = "#B0B0B0", Other = "#E5E5E5")
fig_metrics <- c(Retention_index = "Monocyte retention index", Late_program = "Late lipid/iron programme",
                 adj_TFEC = "TFEC", adj_Iron_handling = "Iron handling")
thf <- theme_classic(base_size = 9) +
  theme(axis.text = element_text(color = "black"), strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9), plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        plot.tag = element_text(face = "bold", size = 13))
star_of <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
get_eff <- function(tp, m) {
  mm <- meta[meta$type == tp & meta$metric == m & is.finite(meta$k) & meta$k > 1, ]
  if (nrow(mm)) return(data.frame(g = mm$g, lo = mm$lo, hi = mm$hi, p = mm$p))
  r <- res[res$type == tp & res$metric == m, ]
  if (!nrow(r)) return(NULL)
  data.frame(g = r$g, lo = r$g - 1.96 * sqrt(r$v), hi = r$g + 1.96 * sqrt(r$v), p = r$p)
}

## a 分期汇总
stg <- c(AAA_vs_Control = "Dilated AAA vs Control",
         Large_vs_Small = "Large vs Small", Ruptured_vs_Large = "Ruptured vs Large")
sa <- bind_rows(lapply(names(stg), function(tp) bind_rows(lapply(names(fig_metrics), function(m) {
  e <- get_eff(tp, m); if (is.null(e)) return(NULL); cbind(type = tp, metric = m, e) }))))
sa$stage <- factor(stg[sa$type], rev(stg)); sa$metric <- factor(fig_metrics[sa$metric], fig_metrics)
sa$star <- star_of(sa$p)
pa <- ggplot(sa, aes(g, stage, color = g > 0)) + geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y", linewidth = 0.7) + geom_point(size = 3) +
  geom_text(aes(x = hi, label = star), hjust = -0.25, vjust = 0.75, size = 5, color = "black") +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") + scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  coord_cartesian(clip = "off") + scale_color_manual(values = c(`TRUE` = col_up, `FALSE` = col_dn), guide = "none") +
  thf + labs(x = "Hedges' g (95% CI)", y = NULL) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13),
        strip.text = element_text(size = 13))

## b / c 单指标森林图（各队列 + 合并）
forest_one <- function(tp, m, title) {
  r <- res[res$type == tp & res$metric == m, ]
  r <- data.frame(what = sprintf("%s (n = %d / %d)", r$dataset, r$n_case, r$n_ctrl), g = r$g,
                  lo = r$g - 1.96 * sqrt(r$v), hi = r$g + 1.96 * sqrt(r$v), pooled = FALSE)
  e <- get_eff(tp, m); mm <- meta[meta$type == tp & meta$metric == m, ]
  r <- rbind(r, data.frame(what = "Pooled (random effects)", g = e$g, lo = e$lo, hi = e$hi, pooled = TRUE))
  r$what <- factor(r$what, rev(r$what))
  sub <- sprintf("g = %.2f, P = %s, I\u00b2 = %.0f%%", e$g, formatC(e$p, format = "g", digits = 2), mm$I2[1] * 100)
  ggplot(r, aes(g, what)) + geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_errorbar(aes(xmin = lo, xmax = hi, color = pooled), width = 0.2, orientation = "y") +
    geom_point(aes(shape = pooled, color = pooled, size = pooled)) +
    scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = col_up), guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 2, `TRUE` = 3.5), guide = "none") +
    thf + labs(x = "Hedges' g (95% CI)", y = NULL, title = title, subtitle = sub) +
    theme(plot.subtitle = element_text(size = 8, hjust = 0.5))
}
pb  <- forest_one("AAA_vs_Control", "Retention_index", "Monocyte retention index\nDilated AAA vs Control")
pc1 <- forest_one("Large_vs_Small", "Retention_index", "Monocyte retention index\nLarge vs Small AAA")
pc2 <- forest_one("Large_vs_Small", "Late_program", "Late lipid/iron programme\nLarge vs Small AAA")

## d GSE98278 分期
d98 <- sc_all[sc_all$dataset == "GSE98278", ]
d98$stage <- factor(ifelse(d98$size == "ruptured", "Ruptured", ifelse(d98$size == "large", "Large", "Intermediate")), names(stage_cols))
ld <- bind_rows(lapply(names(fig_metrics), function(m) data.frame(metric = m, stage = d98$stage, value = d98[[m]])))
pv <- ld %>% group_by(metric) %>% summarise(p = suppressWarnings(wilcox.test(value[stage == "Ruptured"], value[stage == "Large"])$p.value),
                                            ymax = max(value, na.rm = TRUE), .groups = "drop")
pv$lab <- sprintf("P = %.3f", pv$p); pv$metric <- factor(fig_metrics[pv$metric], fig_metrics)
ld$metric <- factor(fig_metrics[ld$metric], fig_metrics)
pd <- ggplot(ld, aes(stage, value, fill = stage)) + geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 0.9) + facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  geom_segment(data = pv, aes(x = 2, xend = 3, y = ymax * 1.08 + 0.05, yend = ymax * 1.08 + 0.05), inherit.aes = FALSE) +
  geom_text(data = pv, aes(x = 2.5, y = ymax * 1.08 + 0.05, label = lab), vjust = -0.5, size = 2.7, inherit.aes = FALSE) +
  scale_fill_manual(values = stage_cols, guide = "none") + thf +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + labs(x = NULL, y = "Score", title = "GSE98278")

## e TFEC 细胞来源
pe <- NULL
if (file.exists(src_csv)) {
  src <- read.csv(src_csv)
  src$lineage <- ifelse(src$celltype %in% names(lineage_map), lineage_map[src$celltype], "Other")
  sl <- src %>% group_by(dataset, lineage) %>% summarise(pct = sum(pct), .groups = "drop")
  sl$lineage <- factor(sl$lineage, names(lin_cols))
  pe <- ggplot(sl, aes(dataset, pct, fill = lineage)) + geom_col(width = 0.65, color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(pct >= 8, sprintf("%.0f%%", pct), "")), position = position_stack(vjust = 0.5), size = 2.7) +
    scale_fill_manual(values = lin_cols, name = NULL) + scale_y_continuous(expand = c(0, 0)) + thf +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = NULL, y = "% of tissue TFEC transcripts", title = "Cellular source of TFEC")
}

## f 小鼠 Tfec 随剂量
pf <- NULL
if (file.exists(dose_csv)) {
  dt <- read.csv(dose_csv); dt <- dt[dt$n >= 30, ]
  dt$dose <- factor(as.character(dt$dose), names(dose_cols)); dt$state <- factor(dt$state, names(mstate_cols))
  pf <- ggplot(dt, aes(dose, mean, color = state, group = state)) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.12) + geom_line(linewidth = 0.9) + geom_point(size = 2) +
    scale_color_manual(values = mstate_cols, name = NULL) + thf + theme(legend.position = "bottom") +
    labs(x = "LiCl (mg/kg)", y = "Tfec expression (mean \u00b1 SE)", title = "Mouse macrophages")
}

## 主图
row1 <- pa + labs(tag = "a")
row2 <- (pb + labs(tag = "b") | pc1 + labs(tag = "c") | pc2) + plot_layout(widths = c(1, 1, 1))
row3 <- pd + labs(tag = "d")
r4 <- Filter(Negate(is.null), list(if (!is.null(pe)) pe + labs(tag = "e"), if (!is.null(pf)) pf + labs(tag = "f")))
rows <- list(row1, row2, row3); hts <- c(0.8, 1, 1)
if (length(r4)) { rows <- c(rows, list(wrap_plots(r4, nrow = 1))); hts <- c(hts, 1.1) }
fig <- wrap_plots(rows, ncol = 1, heights = hts)
pdf_save(fig, "Figure_clinical_main.pdf", 13, 15)
## 各面板单独保存
pdf_save(pa, "Fig_a_stage_summary.pdf", 13, 3.2); pdf_save(pb, "Fig_b_AAA_vs_Control_retention.pdf", 4.5, 3)
pdf_save(pc1, "Fig_c1_Large_vs_Small_retention.pdf", 4.5, 3); pdf_save(pc2, "Fig_c2_Large_vs_Small_late_programme.pdf", 4.5, 3); pdf_save(pd, "Fig_d_GSE98278.pdf", 10, 3.4)
if (!is.null(pe)) pdf_save(pe, "Fig_e_TFEC_source.pdf", 4, 3.6)
if (!is.null(pf)) pdf_save(pf, "Fig_f_mouse_Tfec_dose.pdf", 4, 3.6)

## 补充图：TFEC 敏感性
ss <- sens_sum[sens_sum$type %in% c("AAA_vs_Control", "Neck_vs_Control", "Ruptured_vs_Intact", "Ruptured_vs_Large"), ]
ss$type <- factor(gsub("_", " ", ss$type), rev(gsub("_", " ", c("Neck_vs_Control", "AAA_vs_Control", "Ruptured_vs_Intact", "Ruptured_vs_Large"))))
ps1 <- ggplot(ss, aes(g, type, color = version)) + geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y", position = position_dodge(width = 0.6)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("grey50", "#4292C6", "#B2182B"), name = NULL) + thf + theme(legend.position = "bottom") +
  labs(x = "Hedges' g (95% CI)", y = NULL, title = "TFEC: adjustment for myeloid and B-cell content")
pdf_save(ps1, "Supp_TFEC_sensitivity.pdf", 7, 4)
## 补充图：逐一剔除
lo_df <- loo; lo_df$metric <- factor(primary_lab[lo_df$metric], primary_lab)
full <- meta[meta$type == "AAA_vs_Control" & meta$metric %in% primary, ]
full$metric <- factor(primary_lab[full$metric], primary_lab)
ps2 <- ggplot(lo_df, aes(g, paste("without", dropped))) + geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_vline(data = full, aes(xintercept = g), color = col_up, linewidth = 0.5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y") + geom_point(size = 2) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x") + thf +
  labs(x = "Pooled Hedges' g (95% CI); red line = all 4 cohorts", y = NULL, title = "Leave-one-out: Dilated AAA vs Control")
pdf_save(ps2, "Supp_leave_one_out.pdf", 12, 3)

message("完成，结果在：", out_dir)
