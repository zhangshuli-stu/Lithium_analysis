## =============================================================================
## Monocle3 analysis of the monocyte-macrophage trajectory
##
## Input  mono_mac_5state.rds, carrying sling_pt1 from slingshot_monomac.R
##        (only needed for section 7, the agreement with Slingshot)
##
## Output results_monocle/
##          fig_graph.pdf              principal graph on the UMAP, by state
##          fig_pseudotime.pdf         the same cells by pseudotime
##          fig_modules.pdf            module heatmap + GO terms on bars
##          fig_mhc2.pdf               the MHC II window, pooled and per library
##          fig_profiles.pdf           module mean curves across settings
##          fig_agreement_*.pdf        Slingshot vs Monocle3
##          pseudotime_quartiles.csv, sensitivity.csv, graph_test.csv,
##          modules_k*.csv, enrichment_k*.csv, module_profiles.csv,
##          mhc2_stats.csv, method_agreement.csv, cells.csv
##        plus monocle_pt / monocle_partition on mono_mac
##
## Usage
##   mono_mac <- readRDS("mono_mac_5state.rds")
##   source("monocle3_monomac.R")
##   saveRDS(mono_mac, "mono_mac_5state.rds")   # to keep the new columns
## =============================================================================

library(Seurat)
library(monocle3)
library(ggplot2)
library(patchwork)


# ---- settings ---------------------------------------------------------------

OUT <- "results_monocle"

## trajectory
DROP_FIT   <- c("Resident Mac")        # excluded from the graph entirely
DROP_PT    <- c("Osteoclast-like Mac") # in the graph, out of pseudotime work
ROOT_STATE <- "Classical Mono"         # ASSUMED, from marker identity
NPCS       <- 30
MIN_CELLS  <- 20
SEED       <- 42

## learn_graph grid; the first entry produces the figures
GRID <- list(list(minimal_branch_len = 10),
             list(minimal_branch_len = 15),
             list(minimal_branch_len = 10, ncenter = 200),
             list(minimal_branch_len = 20, ncenter = 300))

## gene modules
MORANS_CUT <- 0.20
K_VALUES   <- c(2, 3)
MAIN_K     <- 3
N_BINS     <- 50
MIN_BIN    <- 10
MIN_DETECT <- 10
TOP_TERMS  <- 8
SIMPLIFY   <- 0.7

## the MHC II panel: the module's core plus the two genes that track it but
## fall below MORANS_CUT (H2-DMa 0.18, Ciita 0.13)
MHC2 <- c("H2-Aa", "H2-Eb1", "H2-Ab1", "Cd74", "H2-DMb1",
          "Tmem176a", "Tmem176b")
MHC1 <- c("H2-K1", "H2-D1", "B2m")     # negative control
MHC2_EXTRA <- c("H2-DMa", "Ciita")     # reported, not clustered

## agreement with Slingshot
SLING_PT_MAX <- 34   # beyond this the Slingshot curve folds back

## Lymphocyte annotations in a myeloid-only dataset come from gene overlap in
## the ontology, not from lymphocyte activity. Flagged, never deleted.
SUSPECT <- c("T cell", "B cell", "lymphocyte", "Natural killer",
             "Th1 and Th2", "Th17", "Primary immunodeficiency", "Allograft",
             "Graft-versus-host", "Type I diabetes", "Asthma",
             "Systemic lupus", "Rheumatoid arthritis")

## appearance
STATE_LEVELS <- c("Classical Mono", "Transitional Mac", "Resident Mac",
                  "Remodeling Mac", "Osteoclast-like Mac")
COLS <- c("Classical Mono"      = "#D98C4A",
          "Transitional Mac"    = "#AFD8E2",
          "Resident Mac"        = "#7ED0F0",
          "Remodeling Mac"      = "#5A61A8",
          "Osteoclast-like Mac" = "#111455")
LIB_COLS  <- c("50" = "#FFD966", "100" = "#6FA8DC",
               "150" = "#8E7CC3", "200" = "#E06666")
MOD_COLS  <- c("#8E9CD4", "#F2C48D", "#8FBF9F", "#D98C8C")
HEAT_COLS <- c("#2B4C9B", "#6FA8D6", "#EAF0F4", "#F6E8A0", "#E8A03C", "#C0392B")
PT_SIZE <- 1.2
FONT    <- "sans"

dir.create(OUT, showWarnings = FALSE)
set.seed(SEED)

void <- theme_void(base_family = FONT, base_size = 10) +
  theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))


# ---- 1. trajectory ----------------------------------------------------------

obj <- mono_mac[, !(mono_mac$myeloid_subtype_short %in% DROP_FIT)]
n    <- table(droplevels(factor(obj$myeloid_subtype_short)))
keep <- names(n)[n >= MIN_CELLS]
if (length(setdiff(names(n), keep))) {
  cat("below", MIN_CELLS, "cells, dropped:",
      paste(setdiff(names(n), keep), collapse = ", "), "\n")
}
obj   <- obj[, obj$myeloid_subtype_short %in% keep]
state <- factor(as.character(obj$myeloid_subtype_short),
                levels = STATE_LEVELS[STATE_LEVELS %in% keep])
lib   <- as.character(obj$sample)

cat(ncol(obj), "cells,", nlevels(state), "states\n"); print(table(state))

cds <- new_cell_data_set(
  GetAssayData(obj, assay = "RNA", layer = "counts"),
  cell_metadata = data.frame(state = as.character(state), library = lib,
                             row.names = colnames(obj)),
  gene_metadata = data.frame(gene_short_name = rownames(obj),
                             row.names = rownames(obj)))
cds <- preprocess_cds(cds, num_dim = NPCS)

# The embedding comes from the Seurat object, so the graph is fitted in the
# same space as the Slingshot run and every figure lines up.
reducedDims(cds)$PCA  <- Embeddings(obj, "harmony")[, 1:NPCS]
reducedDims(cds)$UMAP <- Embeddings(obj, "umap")
cds <- cluster_cells(cds, reduction_method = "UMAP", random_seed = SEED)

cat("\npartitions x state:\n")
print(table(partitions(cds), colData(cds)$state))

# learn_graph overwrites any graph on the object, so every fit starts from
# this untouched copy; otherwise the sensitivity rows describe different
# starting conditions.
cds_clean <- cds
fit <- function(ctrl, base = cds_clean) {
  learn_graph(base, use_partition = TRUE, learn_graph_control = ctrl,
              verbose = FALSE)
}
cds <- fit(GRID[[1]])

root_vertex <- function(cds, root_state) {
  v <- which(colData(cds)$state == root_state)
  closest <- as.matrix(
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)[colnames(cds), ]
  # Confine the root to the main partition: a root in a small side partition
  # leaves the rest of the graph unreachable.
  main <- names(which.max(table(partitions(cds))))
  v    <- v[partitions(cds)[colnames(cds)[v]] == main]
  paste0("Y_", names(which.max(table(closest[v]))))
}
cds <- order_cells(cds, root_pr_nodes = root_vertex(cds, ROOT_STATE))

pt  <- pseudotime(cds)
fin <- is.finite(pt)

quart <- t(sapply(split(pt[fin], droplevels(state[fin])),
                  quantile, c(.25, .5, .75)))
quart <- data.frame(state = rownames(quart),
                    n = as.integer(table(droplevels(state[fin]))[rownames(quart)]),
                    round(quart, 1), row.names = NULL, check.names = FALSE)
cat("\n== pseudotime quartiles ==\n"); print(quart, row.names = FALSE)
write.csv(quart, file.path(OUT, "pseudotime_quartiles.csv"), row.names = FALSE)

flat <- quart$state[quart$`75%` - quart$`25%` < 0.1]
if (length(flat)) {
  cat("\ninterquartile range under 0.1:", paste(flat, collapse = ", "),
      "\n  every cell of that state shares one graph vertex; the value is a\n",
      "  landing point, not a pseudotime, and is excluded below\n")
}
if (sum(!fin)) {
  cat("\n", sum(!fin), "cells unreachable from the root:\n")
  print(table(droplevels(state[!fin])))
}

sens <- do.call(rbind, lapply(seq_along(GRID), function(i) {
  pg <- principal_graph(fit(GRID[[i]]))[["UMAP"]]
  data.frame(setting = i,
             minimal_branch_len = GRID[[i]]$minimal_branch_len,
             ncenter = ifelse(is.null(GRID[[i]]$ncenter), NA, GRID[[i]]$ncenter),
             vertices = igraph::vcount(pg), edges = igraph::ecount(pg),
             branch_points = sum(igraph::degree(pg) > 2),
             leaves = sum(igraph::degree(pg) == 1))
}))
cat("\n== graph structure across learn_graph settings ==\n")
print(sens, row.names = FALSE)
write.csv(sens, file.path(OUT, "sensitivity.csv"), row.names = FALSE)


# ---- 2. trajectory figures --------------------------------------------------
# plot_cells draws the principal graph as flat segment layers that vanish over
# the dark Remodeling cluster. halo() lifts them out and lays a white copy
# underneath.

halo <- function(p, colour = "grey45", lwd = 0.75, mult = 2.2) {
  is_seg <- vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"),
                   logical(1))
  if (!any(is_seg)) return(p)
  seg <- p$layers[is_seg]; p$layers <- p$layers[!is_seg]
  restyle <- function(col, w) lapply(seg, function(l) {
    l2 <- ggproto(NULL, l)
    l2$aes_params$colour <- col; l2$aes_params$linewidth <- w; l2
  })
  p$layers <- c(p$layers, restyle("white", lwd * mult), restyle(colour, lwd))
  p
}

cell_theme <- theme_void(base_family = FONT, base_size = 11) +
  theme(plot.margin = margin(8, 8, 8, 8),
        legend.text = element_text(family = FONT, size = 9),
        legend.title = element_text(family = FONT, size = 9))

p_graph <- halo(
  plot_cells(cds, color_cells_by = "state", label_cell_groups = FALSE,
             label_groups_by_cluster = FALSE, label_leaves = FALSE,
             label_branch_points = FALSE, graph_label_size = 4,
             cell_size = PT_SIZE) +
    scale_colour_manual(values = COLS, drop = FALSE) +
    coord_equal() + cell_theme +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.6))))

p_pt <- halo(
  plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
             label_leaves = FALSE, label_branch_points = FALSE,
             graph_label_size = 4, cell_size = PT_SIZE) +
    coord_equal() + cell_theme + theme(legend.position = "right"))

ggsave(file.path(OUT, "fig_graph.pdf"), p_graph,
       width = 5.6, height = 5.4, device = cairo_pdf)
ggsave(file.path(OUT, "fig_pseudotime.pdf"), p_pt,
       width = 5.8, height = 5, device = cairo_pdf)

umap <- reducedDims(cds)$UMAP
write.csv(data.frame(cell = colnames(cds), state = as.character(state),
                     library = lib, partition = as.character(partitions(cds)),
                     pseudotime = pt, UMAP1 = umap[, 1], UMAP2 = umap[, 2]),
          file.path(OUT, "cells.csv"), row.names = FALSE)


# ---- 3. graph_test ----------------------------------------------------------
# Moran's I on the principal graph, over the whole transcriptome. Fast because
# it runs on the ~166 graph vertices, not on every cell.

gt <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 8)
gt <- gt[order(-gt$morans_I), ]
write.csv(gt, file.path(OUT, "graph_test.csv"))
cat("\n== top 20 by Moran's I ==\n")
print(utils::head(gt[, c("gene_short_name", "morans_I", "q_value")], 20))
cat("\ngenes with a value:", sum(!is.na(gt$morans_I)),
    "| above", MORANS_CUT, ":", sum(gt$morans_I > MORANS_CUT, na.rm = TRUE), "\n")


# ---- 4. pseudotime bins -----------------------------------------------------

use  <- fin & !(as.character(state) %in% DROP_PT)
pt_u <- pt[use]; st_u <- as.character(state)[use]; lb_u <- lib[use]
norm <- normalized_counts(cds)[, use]
cat("\n", sum(use), "cells for module and MHC II analysis\n")

bin   <- cut(pt_u, breaks = N_BINS, labels = FALSE)
good  <- as.integer(names(table(bin))[table(bin) >= MIN_BIN])
bin_pt <- sapply(good, function(b) mean(pt_u[bin == b]))
cat(length(good), "of", N_BINS, "bins have at least", MIN_BIN, "cells\n")

bin_mean <- function(genes) {
  m <- sapply(good, function(b)
    Matrix::rowMeans(norm[genes, bin == b, drop = FALSE]))
  if (length(genes) == 1)
    m <- matrix(m, nrow = 1, dimnames = list(genes, NULL))
  m
}
scale_rows <- function(m) t(apply(m, 1, function(x) {
  r <- range(x); if (diff(r) == 0) rep(0.5, length(x)) else (x - r[1]) / diff(r)
}))


# ---- 5. gene modules and enrichment -----------------------------------------

counts_all <- GetAssayData(mono_mac, assay = "RNA", layer = "counts")
bg <- intersect(rownames(gt)[!is.na(gt$morans_I)],
                rownames(counts_all)[Matrix::rowSums(counts_all > 0) >= MIN_DETECT])
cat(length(bg), "genes in the enrichment background\n")

to_entrez <- function(x) suppressMessages(
  clusterProfiler::bitr(x, "SYMBOL", "ENTREZID", org.Mm.eg.db::org.Mm.eg.db))
bg_map <- to_entrez(bg)

enrich <- function(genes, src) {
  map <- to_entrez(genes)
  if (nrow(map) < 10) return(NULL)
  e <- try({
    if (src == "KEGG")
      clusterProfiler::enrichKEGG(map$ENTREZID, "mmu", universe = bg_map$ENTREZID,
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    else {
      go <- clusterProfiler::enrichGO(map$ENTREZID, org.Mm.eg.db::org.Mm.eg.db,
                                      ont = "BP", universe = bg_map$ENTREZID,
                                      pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                                      readable = TRUE)
      # GO BP returns one process at several levels of its hierarchy.
      if (!is.null(go) && nrow(as.data.frame(go)) > 1)
        go <- clusterProfiler::simplify(go, cutoff = SIMPLIFY,
                                        by = "p.adjust", select_fun = min)
      go
    }}, silent = TRUE)
  if (inherits(e, "try-error") || is.null(e)) return(NULL)
  d <- as.data.frame(e); if (nrow(d)) d else NULL
}

sig    <- intersect(rownames(gt)[which(gt$morans_I > MORANS_CUT)], rownames(norm))
scaled <- scale_rows(bin_mean(sig))
cat("\n", length(sig), "genes with Moran's I >", MORANS_CUT, "\n")

hc       <- hclust(stats::as.dist(1 - cor(t(scaled))), method = "ward.D2")
profiles <- list(); results <- list()

for (K in K_VALUES) {
  cl <- cutree(hc, k = K)
  # Number modules by where they peak, so module 1 is early for every K.
  peak <- sapply(seq_len(K), function(m) weighted.mean(
    seq_along(bin_pt), colMeans(scaled[cl == m, , drop = FALSE])))
  mod  <- as.integer(setNames(rank(peak, ties.method = "first"),
                              seq_len(K))[as.character(cl)])
  
  # Coherence: mean correlation of each gene with its own module's profile.
  # A module well below 0.5 is a bin of leftover shapes, not a programme.
  coh <- sapply(seq_len(K), function(m) {
    g <- scaled[mod == m, , drop = FALSE]
    if (nrow(g) < 2) NA_real_ else mean(apply(g, 1, function(x) cor(x, colMeans(g))))
  })
  
  mods <- data.frame(gene = rownames(scaled), module = mod,
                     morans_I = gt[rownames(scaled), "morans_I"],
                     row.names = NULL)
  mods <- mods[order(mods$module, -mods$morans_I), ]
  write.csv(mods, file.path(OUT, sprintf("modules_k%d.csv", K)), row.names = FALSE)
  
  cat("\n-- k =", K, "--\n")
  print(data.frame(module = seq_len(K), n = as.integer(table(mod)),
                   coherence = round(coh, 2),
                   peak_pt = round(bin_pt[pmax(1, round(peak))], 1)),
        row.names = FALSE)
  
  profiles[[as.character(K)]] <- do.call(rbind, lapply(seq_len(K), function(m) {
    g <- scaled[mod == m, , drop = FALSE]
    data.frame(k = K, module = m, pt = bin_pt, mean = colMeans(g),
               n = nrow(g), coherence = coh[m])
  }))
  
  enr <- list()
  for (m in seq_len(K)) for (src in c("GOBP", "KEGG")) {
    d <- enrich(mods$gene[mods$module == m], src)
    if (is.null(d)) next
    d$module <- m; d$source <- src
    d$suspect <- grepl(paste(SUSPECT, collapse = "|"), d$Description,
                       ignore.case = TRUE)
    enr[[length(enr) + 1]] <- d[, c("module", "source", "ID", "Description",
                                    "GeneRatio", "geneID", "pvalue",
                                    "p.adjust", "Count", "suspect")]
  }
  if (!length(enr)) next
  enr <- do.call(rbind, enr); enr <- enr[order(enr$module, enr$pvalue), ]
  write.csv(enr, file.path(OUT, sprintf("enrichment_k%d.csv", K)), row.names = FALSE)
  results[[as.character(K)]] <- list(K = K, mods = mods, enr = enr, coh = coh)
  
  for (m in seq_len(K)) {
    sub <- utils::head(enr[enr$source == "GOBP" & enr$module == m, ], TOP_TERMS)
    if (!nrow(sub)) next
    cat("  module", m, "- GO BP:\n")
    for (i in seq_len(nrow(sub)))
      cat(sprintf("    %-52s p.adj=%.1e n=%d%s\n",
                  substr(sub$Description[i], 1, 52), sub$p.adjust[i],
                  sub$Count[i], ifelse(sub$suspect[i], "  [lymphocyte]", "")))
  }
}

prof <- do.call(rbind, profiles)
write.csv(prof, file.path(OUT, "module_profiles.csv"), row.names = FALSE)

prof$panel <- paste("k =", prof$k)
ggsave(file.path(OUT, "fig_profiles.pdf"),
       ggplot(prof, aes(pt, mean, colour = factor(module))) +
         geom_line(linewidth = 0.8) +
         scale_colour_manual(values = MOD_COLS, name = "Module") +
         facet_wrap(~ panel, nrow = 1) +
         labs(x = "Pseudotime", y = "Mean scaled expression") +
         theme_classic(base_family = FONT, base_size = 10) +
         theme(strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 9)),
       width = 3.3 * length(K_VALUES), height = 3, device = cairo_pdf)


# ---- 6. module figure -------------------------------------------------------
# Heatmap, a module colour strip, and GO terms written onto the bars. Three
# stitched panels rather than one faceted plot, because the bar axis has to
# sit on top and the labels have to land on the bars.

main <- results[[as.character(MAIN_K)]]
K    <- main$K
mods <- main$mods

ord <- order(mods$module,
             apply(scaled[mods$gene, ], 1,
                   function(x) weighted.mean(seq_along(x), x)))
genes_ord <- mods$gene[ord]

hm <- expand.grid(gene = genes_ord, bin = seq_along(bin_pt),
                  stringsAsFactors = FALSE)
hm$value  <- scaled[cbind(match(hm$gene, rownames(scaled)), hm$bin)]
hm$module <- factor(mods$module[match(hm$gene, mods$gene)])
hm$gene   <- factor(hm$gene, levels = rev(genes_ord))

top <- do.call(rbind, lapply(seq_len(K), function(m) utils::head(
  main$enr[main$enr$source == "GOBP" & main$enr$module == m, ], TOP_TERMS)))
top$module <- factor(top$module, levels = seq_len(K))
top$logp   <- -log10(top$p.adjust)
top$Description <- sub(" of exogenous peptide antigen via", " via", top$Description)
top$Description <- sub(" with MHC class II protein complex", " with MHC II",
                       top$Description)

# y counted downwards with a gap between modules, so heatmap blocks and bar
# blocks line up row for row.
top$y <- NA_real_; ypos <- 0
for (m in levels(top$module)) {
  i <- which(top$module == m)
  top$y[i] <- ypos + seq_along(i); ypos <- max(top$y[i]) + 0.8
}
top$y <- -top$y
ylim  <- range(top$y) + c(-0.9, 0.9)

blocks <- do.call(rbind, lapply(levels(top$module), function(m) {
  d <- top[top$module == m, ]
  data.frame(module = factor(m, levels = levels(top$module)),
             ymin = min(d$y) - 0.5, ymax = max(d$y) + 0.5)
}))

strip <- data.frame(gene = factor(genes_ord, levels = rev(genes_ord)),
                    module = factor(mods$module[match(genes_ord, mods$gene)]))
p_strip <- ggplot(strip, aes(1, gene, fill = module)) +
  geom_raster() + scale_fill_manual(values = MOD_COLS[seq_len(K)]) +
  scale_y_discrete(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  void

ticks <- round(seq(1, length(bin_pt), length.out = 4))
# geom_tile, not geom_raster: raster is a bitmap and 164 x 48 cells in a
# narrow panel blurs. Tiles stay vector at the cost of a larger file.
p_hm <- ggplot(hm, aes(bin, gene, fill = value)) +
  geom_tile() +
  scale_x_continuous(breaks = ticks, labels = round(bin_pt[ticks], 1),
                     expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(colours = HEAT_COLS, name = NULL, breaks = c(0, 0.5, 1),
                       guide = guide_colourbar(direction = "horizontal",
                                               barwidth = unit(38, "mm"),
                                               barheight = unit(3, "mm"),
                                               ticks = FALSE)) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Trajectory", x = NULL, y = NULL) +
  theme_void(base_family = FONT, base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, margin = margin(b = 3)),
        strip.text = element_blank(),
        axis.text.x = element_text(size = 8.5, colour = "#444444",
                                   margin = margin(t = 2)),
        legend.position = "bottom", legend.text = element_text(size = 7),
        legend.margin = margin(t = 1), plot.margin = margin(2, 2, 2, 2))

xmax <- max(top$logp) * 1.06
p_bar <- ggplot(top) +
  geom_rect(data = blocks, aes(xmin = 0, xmax = xmax, ymin = ymin, ymax = ymax,
                               fill = module), alpha = 0.16) +
  geom_col(aes(logp, y, fill = module), width = 0.62, orientation = "y") +
  geom_text(aes(xmax * 0.012, y, label = Description), hjust = 0, size = 3.6,
            family = FONT, colour = "#1A1A1A") +
  # Asterisks mark lymphocyte-specific terms: MHC II's own annotation includes
  # T cell activation, since that is what it does, but there are no T cells here.
  geom_text(data = top[top$suspect, ], aes(xmax * 0.012, y), label = "*",
            hjust = 1.4, size = 5, colour = "#B03A2E", family = FONT) +
  scale_fill_manual(values = MOD_COLS[seq_len(K)]) +
  scale_x_continuous(position = "top", limits = c(0, xmax), expand = c(0, 0),
                     breaks = scales::pretty_breaks(4)) +
  scale_y_continuous(limits = ylim, expand = c(0, 0)) +
  labs(x = expression(-log[10] ~ italic(P)[adj])) +
  theme_void(base_family = FONT, base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9, colour = "#333333"),
        axis.title.x.top = element_text(size = 10, hjust = 0, margin = margin(b = 2)),
        axis.ticks.x.top = element_line(colour = "#333333", linewidth = 0.3),
        axis.ticks.length.x.top = unit(1.2, "mm"),
        plot.margin = margin(2, 4, 2, 2))

leg <- data.frame(module = factor(seq_len(K)), lab = paste("Cluster", seq_len(K)))
p_leg <- ggplot(leg, aes(seq_len(K), 1, fill = module)) +
  geom_point(shape = 22, size = 4, colour = NA) +
  geom_text(aes(label = lab), nudge_x = 0.13, hjust = 0, size = 3.8, family = FONT) +
  scale_fill_manual(values = MOD_COLS[seq_len(K)]) +
  scale_x_continuous(limits = c(0.6, K + 1.2), expand = c(0, 0)) + void

ggsave(file.path(OUT, "fig_modules.pdf"),
       p_leg / (p_strip + p_hm + p_bar +
                  plot_layout(widths = c(0.035, 0.75, 1.15), nrow = 1)) +
         plot_layout(heights = c(0.045, 1)),
       width = 8.6, height = 0.32 * nrow(top) + 1.6, device = cairo_pdf)


# ---- 7. the MHC II window ---------------------------------------------------
# Raw expression, not scaled: the size of the change is the point, and
# row-scaling would map a 4-fold and a 1.2-fold rise onto the same 0-1 range.

mhc_all <- intersect(c(MHC2, MHC2_EXTRA, MHC1), rownames(norm))
binned  <- bin_mean(mhc_all)
win     <- range(bin_pt[which(scale_rows(binned)["Cd74", ] > 0.6)])
in_win  <- pt_u >= win[1] & pt_u <= win[2]

fold <- sapply(mhc_all, function(g)
  mean(norm[g, in_win]) / mean(norm[g, !in_win]))
detect <- sapply(mhc_all, function(g)
  c(in_window = mean(norm[g, in_win] > 0), outside = mean(norm[g, !in_win] > 0)))

mhc_stats <- data.frame(gene = mhc_all,
                        set = ifelse(mhc_all %in% MHC1, "MHC I control",
                                     ifelse(mhc_all %in% MHC2_EXTRA,
                                            "below Moran cutoff", "module")),
                        fold = round(fold, 2),
                        cor_with_Cd74 = round(
                          apply(binned, 1, function(x) cor(x, binned["Cd74", ])), 2),
                        detected_in = round(detect["in_window", ], 3),
                        detected_out = round(detect["outside", ], 3),
                        row.names = NULL)
cat(sprintf("\n== MHC II window: pseudotime %.1f to %.1f, %d cells ==\n",
            win[1], win[2], sum(in_win)))
print(table(st_u[in_win]))
cat("share of each library's cells inside the window:\n")
print(round(table(lb_u[in_win]) / table(lb_u), 3))
print(mhc_stats, row.names = FALSE)
write.csv(mhc_stats, file.path(OUT, "mhc2_stats.csv"), row.names = FALSE)

band <- annotate("rect", xmin = win[1], xmax = win[2], ymin = -Inf, ymax = Inf,
                 fill = "#F0C060", alpha = 0.18)
mhc_df <- data.frame(
  gene  = factor(rep(MHC2, times = length(bin_pt)), levels = MHC2),
  pt    = rep(bin_pt, each = length(MHC2)),
  value = as.vector(binned[MHC2, ]))

# Per library, to show the peak is not carried by one of them. Half as many
# bins, since a single library has far fewer cells per bin than the pool.
mhc_lib <- do.call(rbind, lapply(sort(unique(lb_u)), function(L) {
  s  <- lb_u == L
  b2 <- cut(pt_u[s], breaks = N_BINS / 2, labels = FALSE)
  ok <- as.integer(names(table(b2))[table(b2) >= MIN_BIN])
  do.call(rbind, lapply(MHC2, function(g) data.frame(
    gene = g, library = L,
    pt    = sapply(ok, function(b) mean(pt_u[s][b2 == b])),
    value = sapply(ok, function(b) mean(norm[g, s][b2 == b])))))
}))
mhc_lib$gene <- factor(mhc_lib$gene, levels = MHC2)

line_theme <- theme_classic(base_family = FONT, base_size = 10) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9))

ggsave(file.path(OUT, "fig_mhc2.pdf"),
       (ggplot(mhc_df, aes(pt, value)) + band +
          geom_line(linewidth = 0.7, colour = "#2B2B2B") +
          facet_wrap(~ gene, scales = "free_y", nrow = 2) +
          labs(x = "Pseudotime", y = "Mean expression") + line_theme) /
         (ggplot(mhc_lib, aes(pt, value, colour = library)) + band +
            geom_line(linewidth = 0.8) +
            scale_colour_manual(values = LIB_COLS, name = "mg/kg") +
            facet_wrap(~ gene, scales = "free_y", nrow = 2) +
            labs(x = "Pseudotime", y = "Mean expression") + line_theme +
            theme(legend.position = "bottom")) +
         plot_layout(heights = c(1, 1.2)),
       width = 9.5, height = 8, device = cairo_pdf)


# ---- 8. agreement with Slingshot --------------------------------------------

if ("sling_pt1" %in% colnames(mono_mac@meta.data)) {
  s_pt <- mono_mac$sling_pt1[names(pt)]
  lin  <- setdiff(levels(state), DROP_PT)
  ag   <- data.frame(sling = s_pt, mono = pt,
                     state = factor(as.character(state), levels = lin))
  ag   <- ag[!is.na(ag$state) & is.finite(ag$sling) & is.finite(ag$mono), ]
  ag$kept <- ag$sling <= SLING_PT_MAX
  
  rho_all  <- cor(ag$sling, ag$mono, method = "spearman")
  rho_kept <- cor(ag$sling[ag$kept], ag$mono[ag$kept], method = "spearman")
  rho_st   <- sapply(lin, function(s) {
    i <- ag$kept & ag$state == s
    if (sum(i) < 50) NA_real_ else cor(ag$sling[i], ag$mono[i], method = "spearman")
  })
  cat(sprintf("\n== agreement: all %d cells rho = %.3f | <= %d, %d cells rho = %.3f ==\n",
              nrow(ag), rho_all, SLING_PT_MAX, sum(ag$kept), rho_kept))
  print(round(rho_st, 3))
  
  write.csv(data.frame(
    comparison = c("lineage states, both methods", 
                   sprintf("slingshot pseudotime <= %d", SLING_PT_MAX), lin),
    n = c(nrow(ag), sum(ag$kept), as.integer(table(ag$state[ag$kept])[lin])),
    spearman = round(c(rho_all, rho_kept, rho_st), 3)),
    file.path(OUT, "method_agreement.csv"), row.names = FALSE)
  
  # A scatter would show horizontal banding: Monocle3 projects cells onto ~166
  # vertices, so its pseudotime takes few distinct values. Hex bins instead,
  # with a binned-median trend rather than a loess, which bends down at the
  # boundary.
  kept <- ag[ag$kept, ]
  grp  <- findInterval(kept$sling, seq(0, SLING_PT_MAX, length.out = 31),
                       rightmost.closed = TRUE)
  trend <- do.call(rbind, lapply(split(seq_len(nrow(kept)), grp), function(i) {
    if (length(i) < 20) return(NULL)
    d <- kept[i, ]
    data.frame(sling = median(d$sling), med = median(d$mono),
               lo = as.numeric(quantile(d$mono, .25)),
               hi = as.numeric(quantile(d$mono, .75)))
  }))
  
  lab <- sprintf("rho == %.2f ~ (n == '%s')", rho_kept,
                 format(sum(ag$kept), big.mark = ","))
  ggsave(file.path(OUT, "fig_agreement_density.pdf"),
         ggplot(kept, aes(sling, mono)) +
           geom_hex(bins = 45) +
           scale_fill_gradientn(colours = c("#EAF0F4", "#9EC4DE", "#4C7FB8",
                                            "#2B4C9B", "#16265A"),
                                trans = "log10", name = "Cells",
                                labels = scales::label_number(accuracy = 1)) +
           geom_ribbon(data = trend, aes(sling, ymin = lo, ymax = hi),
                       inherit.aes = FALSE, fill = "#4A4A4A", alpha = 0.14) +
           geom_line(data = trend, aes(sling, med), inherit.aes = FALSE,
                     colour = "#2B2B2B", linewidth = 0.8) +
           annotate("text", x = 0, y = max(kept$mono), hjust = 0, vjust = 1,
                    label = lab, parse = TRUE, size = 3.1, family = FONT) +
           labs(x = "Slingshot pseudotime", y = "Monocle3 pseudotime") +
           theme_classic(base_family = FONT, base_size = 10) +
           theme(legend.position = c(0.97, 0.03),
                 legend.justification = c(1, 0),
                 legend.background = element_blank(),
                 legend.key.width = unit(3, "mm"),
                 legend.key.height = unit(3.5, "mm"),
                 legend.title = element_text(size = 8),
                 legend.text = element_text(size = 7)),
         width = 4.8, height = 4.4, device = cairo_pdf)
  
  kept$bin <- cut(kept$sling, breaks = 28, labels = FALSE)
  agg <- do.call(rbind, lapply(
    split(kept, list(kept$state, kept$bin), drop = TRUE), function(d) {
      if (nrow(d) < 15) return(NULL)
      data.frame(state = d$state[1], sling = median(d$sling),
                 med = median(d$mono),
                 lo = as.numeric(quantile(d$mono, .25)),
                 hi = as.numeric(quantile(d$mono, .75)))
    }))
  ggsave(file.path(OUT, "fig_agreement_states.pdf"),
         ggplot(agg, aes(sling, med, colour = state, fill = state)) +
           geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
           geom_line(linewidth = 1.1) +
           scale_colour_manual(values = COLS, name = NULL) +
           scale_fill_manual(values = COLS, guide = "none") +
           labs(x = "Slingshot pseudotime", y = "Monocle3 pseudotime") +
           guides(colour = guide_legend(override.aes = list(linewidth = 1.6))) +
           theme_classic(base_family = FONT, base_size = 10) +
           theme(legend.position = "bottom"),
         width = 4.8, height = 4.8, device = cairo_pdf)
} else {
  cat("\nsling_pt1 not on the object; section 8 skipped\n")
}


# ---- 9. write back ----------------------------------------------------------
# Matched on barcode. Cells outside the fit get NA rather than a neighbour's
# value. Cells whose state was dropped from pseudotime work are masked.

put <- function(values_named, name) {
  v <- rep(NA_real_, ncol(mono_mac)); names(v) <- colnames(mono_mac)
  hit <- intersect(colnames(mono_mac), names(values_named))
  v[hit] <- values_named[hit]
  mono_mac[[name]] <- v
  mono_mac
}
pt_masked <- pt
pt_masked[as.character(state) %in% DROP_PT | !fin] <- NA_real_
mono_mac <- put(pt_masked, "monocle_pt")
mono_mac <- put(setNames(as.numeric(partitions(cds)), colnames(cds)),
                "monocle_partition")

mono_mac@misc$monocle3 <- list(
  dropped_from_fit = DROP_FIT,
  dropped_from_pseudotime = DROP_PT,
  root = paste(ROOT_STATE, "- ASSUMED from Ccr2/Ly6c2 identity"),
  reduction = sprintf("harmony (%d dims), UMAP from the Seurat object", NPCS),
  morans_cutoff = MORANS_CUT, modules_k = MAIN_K,
  clustering = "Ward linkage on 1 - correlation of binned profiles",
  cells_fitted = ncol(cds), date = as.character(Sys.Date()))

mono_mac <<- mono_mac
cat("\nwrote", OUT, "/\n")
cat("added monocle_pt, monocle_partition to mono_mac\n")
cat("save with saveRDS(mono_mac, \"mono_mac_5state.rds\") to keep them\n")
