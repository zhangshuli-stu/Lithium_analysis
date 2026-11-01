############################################################################
## 转录因子筛选的完整证据链
##   图 1  筛选示意：289 个候选转录因子的人鼠轨迹趋势散点，37 个一致的着色，6 个 key TFs 标名
##   图 2  系统检验：37 个人鼠一致的转录因子全部在 bulk 队列中检验（37 个之间做 FDR 校正）
##         2a 破裂 vs 大瘤体排序图   2b 火山图   2c 四类临床比较的效应量热图
## 全部离线运行；依赖 05_transition 的 P3 结果与 rawdata 中的 bulk 数据
############################################################################
tf_csv   <- "/home/data/t080553/Li_AAA/05_transition/P3_TF_expression_conservation.csv"
bulk_dir <- "/home/data/t080553/Li_AAA/rawdata"
out_dir  <- "/home/data/t080553/Li_AAA/08_TF_screen"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

key_tfs <- c("MAF", "TFEC", "BHLHE41", "NFE2L1", "EPAS1", "SMAD3")
neck_ids <- c("GSE47472")
pan_myeloid <- c("CD68", "CSF1R", "C1QA", "C1QB", "C1QC", "AIF1", "TYROBP", "FCER1G", "CD14", "MS4A7")
b_genes     <- c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1")

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })
filter <- dplyr::filter; mutate <- dplyr::mutate; summarise <- dplyr::summarise; group_by <- dplyr::group_by
select <- dplyr::select; arrange <- dplyr::arrange; n <- dplyr::n
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
th <- theme_classic(base_size = 11) + theme(axis.text = element_text(color = "black"),
                                            plot.title = element_text(face = "bold", hjust = 0.5),
                                            strip.background = element_blank(), strip.text = element_text(face = "bold"))
pdf_save <- function(p, name, w, h) ggsave(file.path(out_dir, name), p, width = w, height = h, bg = "white")
col_up <- "#B2182B"; col_dn <- "#2166AC"

######################################################################
## 图 1：筛选示意
######################################################################
tf <- read.csv(tf_csv)
tf$cls <- ifelse(tf$conserved, ifelse(tf$trend_mouse > 0, "Conserved, up", "Conserved, down"), "Not conserved")
tf$cls <- factor(tf$cls, c("Conserved, up", "Conserved, down", "Not conserved"))
tf$is_key <- tf$tf %in% key_tfs
cat("候选转录因子：", nrow(tf), "；人鼠一致：", sum(tf$conserved), "\n")
cat("key TFs 是否都在一致列表中：\n"); print(setNames(key_tfs %in% tf$tf[tf$conserved], key_tfs))

p1 <- ggplot(tf, aes(trend_mouse, trend_human)) +
  annotate("rect", xmin = 0.5, xmax = 1.05, ymin = 0.5, ymax = 1.05, fill = col_up, alpha = 0.06) +
  annotate("rect", xmin = -1.05, xmax = -0.5, ymin = -1.05, ymax = -0.5, fill = col_dn, alpha = 0.06) +
  geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
  geom_point(data = subset(tf, cls == "Not conserved"), color = "grey80", size = 1.4) +
  geom_point(data = subset(tf, cls != "Not conserved"), aes(color = cls, size = r_mouse_human), alpha = 0.85) +
  geom_point(data = subset(tf, is_key), shape = 21, size = 4, stroke = 1, color = "black", fill = NA) +
  scale_color_manual(values = c("Conserved, up" = col_up, "Conserved, down" = col_dn), name = NULL) +
  scale_size(range = c(1.5, 4), name = "Mouse-human\ncurve r") +
  coord_fixed(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05)) + th +
  labs(x = "Mouse: trend along trajectory (Spearman \u03c1)", y = "Human: trend along trajectory (Spearman \u03c1)",
       title = sprintf("%d of %d transcription factors conserved", sum(tf$conserved), nrow(tf)))
lab_df <- subset(tf, is_key)
p1 <- p1 + if (has_repel) ggrepel::geom_text_repel(data = lab_df, aes(label = tf), fontface = "bold.italic", size = 4,
                                                   box.padding = 0.6, min.segment.length = 0) else
                                                     geom_text(data = lab_df, aes(label = tf), fontface = "bold.italic", size = 4, vjust = -1.2)
pdf_save(p1, "Fig1_TF_screen_scatter.pdf", 6.5, 5.5)

######################################################################
## 图 2：37 个人鼠一致的转录因子在 bulk 队列中的系统检验
######################################################################
cons_tfs <- tf$tf[tf$conserved]
files <- list.files(bulk_dir, pattern = "^bulk_GSE\\d+\\.rds$", full.names = TRUE)
bulk <- lapply(files, function(f) {
  x <- readRDS(f); ex <- as.matrix(x$exp); storage.mode(ex) <- "numeric"
  if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)
  ex <- ex[rowSums(is.na(ex)) == 0 & apply(ex, 1, var) > 0, , drop = FALSE]
  list(id = sub("^bulk_(GSE\\d+)\\.rds$", "\\1", basename(f)), exp = ex, group = tolower(as.character(x$group)),
       size = if (!is.null(x$size)) tolower(as.character(x$size)) else rep(NA, ncol(ex)))
})
names(bulk) <- sapply(bulk, `[[`, "id")
score_set <- function(ex, g) { g <- intersect(g, rownames(ex)); if (length(g) < 3) return(rep(NA, ncol(ex)))
colMeans(t(scale(t(ex[g, , drop = FALSE])))) }
tf_scores <- lapply(bulk, function(b) {
  pan <- score_set(b$exp, pan_myeloid); bs <- score_set(b$exp, b_genes)
  out <- lapply(cons_tfs, function(t) {
    if (!t %in% rownames(b$exp)) return(NULL)
    z <- as.numeric(scale(b$exp[t, ]))
    data.frame(dataset = b$id, sample = colnames(b$exp), group = b$group, size = b$size, tf = t, raw = z,
               adj = resid(lm(z ~ pan)), adj2 = resid(lm(z ~ pan + bs)))
  })
  bind_rows(out)
})
ts <- bind_rows(tf_scores)

hedges_g <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]; n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(c(g = NA, v = NA, p = NA))
  sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  g <- (1 - 3 / (4 * (n1 + n2) - 9)) * (mean(x) - mean(y)) / sp
  c(g = g, v = (n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2)), p = suppressWarnings(wilcox.test(x, y)$p.value))
}
meta_dl <- function(g, v) {
  ok <- is.finite(g) & is.finite(v); g <- g[ok]; v <- v[ok]; k <- length(g)
  if (k == 1) return(c(g = g, lo = g - 1.96 * sqrt(v), hi = g + 1.96 * sqrt(v), p = NA, k = 1))
  if (k < 1) return(c(g = NA, lo = NA, hi = NA, p = NA, k = 0))
  w <- 1 / v; gf <- sum(w * g) / sum(w); Q <- sum(w * (g - gf)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (v + tau2); gr <- sum(ws * g) / sum(ws); se <- sqrt(1 / sum(ws))
  c(g = gr, lo = gr - 1.96 * se, hi = gr + 1.96 * se, p = 2 * pnorm(-abs(gr / se)), k = k)
}
is_ctrl <- function(x) grepl("control|normal|ctrl|donor", x)
## 每个比较：(数据集, case 逻辑, control 逻辑)
cmp_def <- function(d, type) switch(type,
                                    AAA_vs_Control    = if (!d$dataset[1] %in% neck_ids && any(is_ctrl(d$group))) list(!is_ctrl(d$group), is_ctrl(d$group)),
                                    Large_vs_Small    = if (d$dataset[1] == "GSE57691") list(grepl("large", d$group), grepl("small", d$group)) else
                                      if (d$dataset[1] == "GSE98278") list(d$size == "large", d$size == "intermediate"),
                                    Ruptured_vs_Intact = if (d$dataset[1] == "GSE98278") list(d$size == "ruptured", d$size != "ruptured"),
                                    Ruptured_vs_Large  = if (d$dataset[1] == "GSE98278") list(d$size == "ruptured", d$size == "large"))
types <- c("AAA_vs_Control", "Large_vs_Small", "Ruptured_vs_Intact", "Ruptured_vs_Large")

res <- bind_rows(lapply(types, function(tp) bind_rows(lapply(cons_tfs, function(t) bind_rows(lapply(c("raw", "adj", "adj2"), function(vv) {
  per <- bind_rows(lapply(split(ts[ts$tf == t, ], ts$dataset[ts$tf == t]), function(d) {
    cs <- cmp_def(d, tp); if (is.null(cs)) return(NULL)
    h <- hedges_g(d[[vv]][cs[[1]]], d[[vv]][cs[[2]]]); data.frame(dataset = d$dataset[1], g = h["g"], v = h["v"], p = h["p"])
  }))
  if (!nrow(per)) return(NULL)
  if (nrow(per) == 1) return(data.frame(type = tp, tf = t, version = vv, g = per$g, p = per$p, k = 1))
  m <- meta_dl(per$g, per$v); data.frame(type = tp, tf = t, version = vv, g = m["g"], p = m["p"], k = m["k"])
}))))))
res <- res %>% group_by(type, version) %>% mutate(fdr = p.adjust(p, "BH"), rank = rank(p, ties.method = "min")) %>% ungroup() %>% as.data.frame()
write.csv(res, file.path(out_dir, "TF37_bulk_tests.csv"), row.names = FALSE)

cat("\n== TFEC 在 37 个人鼠一致转录因子中的排名（按 P 值；巨噬含量校正版本）==\n")
print(res %>% filter(tf == "TFEC") %>% select(type, version, g, p, fdr, rank) %>%
        mutate(g = round(g, 2), p = signif(p, 2), fdr = signif(fdr, 2)), row.names = FALSE)
cat("\n== 破裂 vs 大瘤体：排名前 10（巨噬含量校正）==\n")
print(res %>% filter(type == "Ruptured_vs_Large", version == "adj") %>% arrange(p) %>% head(10) %>%
        select(tf, g, p, fdr, rank) %>% mutate(g = round(g, 2), p = signif(p, 2), fdr = signif(fdr, 2)), row.names = FALSE)

## 2a 排序图：破裂 vs 大瘤体（巨噬含量校正）
ra <- res %>% filter(type == "Ruptured_vs_Large", version == "adj") %>% arrange(g)
ra$tf <- factor(ra$tf, ra$tf)
ra$sig <- ifelse(ra$p < 0.05, "P < 0.05", "n.s.")
p2a <- ggplot(ra, aes(g, tf, fill = interaction(g > 0, sig))) + geom_col(width = 0.75) +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_manual(values = c("TRUE.P < 0.05" = col_up, "FALSE.P < 0.05" = col_dn, "TRUE.n.s." = "#F4B6B6", "FALSE.n.s." = "#B6CDE8"),
                    labels = c("TRUE.P < 0.05" = "Higher, P < 0.05", "FALSE.P < 0.05" = "Lower, P < 0.05",
                               "TRUE.n.s." = "Higher, n.s.", "FALSE.n.s." = "Lower, n.s."), name = NULL) +
  th + theme(axis.text.y = element_text(face = ifelse(levels(ra$tf) %in% key_tfs, "bold.italic", "italic"),
                                        color = ifelse(levels(ra$tf) == "TFEC", col_dn, "black"), size = 8)) +
  labs(x = "Hedges' g (myeloid-adjusted)", y = NULL, title = "Ruptured vs Large AAA (GSE98278)")
pdf_save(p2a, "Fig2a_TF37_rank_ruptured_vs_large.pdf", 5.5, 7)

## 2b 火山图：两个破裂比较
vb <- res %>% filter(type %in% c("Ruptured_vs_Intact", "Ruptured_vs_Large"), version == "adj")
vb$type <- factor(gsub("_", " ", vb$type), c("Ruptured vs Intact", "Ruptured vs Large"))
vb$hl <- vb$tf %in% key_tfs
p2b <- ggplot(vb, aes(g, -log10(p))) + geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey60") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_point(aes(color = ifelse(p < 0.05, ifelse(g > 0, "up", "down"), "ns")), size = 2) +
  scale_color_manual(values = c(up = col_up, down = col_dn, ns = "grey70"), guide = "none") +
  facet_wrap(~ type) + th + labs(x = "Hedges' g (myeloid-adjusted)", y = expression(-log[10]~P))
p2b <- p2b + if (has_repel) ggrepel::geom_text_repel(data = subset(vb, hl | p < 0.05), aes(label = tf), size = 3.2,
                                                     fontface = "italic", min.segment.length = 0) else
                                                       geom_text(data = subset(vb, hl | p < 0.05), aes(label = tf), size = 3.2, fontface = "italic", vjust = -0.8)
pdf_save(p2b, "Fig2b_TF37_volcano_rupture.pdf", 9, 4.2)

## 2c 热图：四类临床比较（巨噬含量校正）
hc <- res %>% filter(version == "adj")
ord <- hc %>% filter(type == "Ruptured_vs_Large") %>% arrange(g) %>% pull(tf)
hc$tf <- factor(hc$tf, rev(ord))
hc$type <- factor(gsub("_", " ", hc$type), gsub("_", " ", types))
hc$star <- ifelse(is.na(hc$p), "", ifelse(hc$p < 0.001, "***", ifelse(hc$p < 0.01, "**", ifelse(hc$p < 0.05, "*", ""))))
p2c <- ggplot(hc, aes(type, tf, fill = g)) + geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = star), size = 3.2, vjust = 0.75) +
  scale_fill_gradient2(low = col_dn, mid = "white", high = col_up, limits = c(-1.5, 1.5), oob = scales::squish, name = "Hedges' g") +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
             axis.text.y = element_text(face = ifelse(rev(ord) %in% key_tfs, "bold.italic", "italic"), size = 8)) +
  labs(x = NULL, y = NULL, title = "Conserved TFs across clinical comparisons")
pdf_save(p2c, "Fig2c_TF37_heatmap_clinical.pdf", 5.2, 7.5)

## 组合
pdf_save((p1 | p2a) / p2b + plot_layout(heights = c(1.4, 1)), "Fig_TF_screen_combined.pdf", 13, 11)
message("完成，结果在：", out_dir)
