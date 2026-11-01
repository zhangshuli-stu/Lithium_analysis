## =============================================================================
## Slingshot trajectory analysis of the monocyte-macrophage compartment
##
## Input :  mono_mac_5state.rds   (16,945 cells, 5 states, QC-filtered)
## Output:  results_slingshot/
##            fig_states.pdf          UMAP by state, with lineage arrows
##            fig_pseudotime.pdf      UMAP coloured by lineage-1 pseudotime
##            fig_combined.pdf        the two side by side
##            cells_<TAG>.csv         per-cell pseudotime and curve weights
##            curves_<TAG>.csv        curve coordinates in UMAP space
##            pseudotime_quartiles_<TAG>.csv
##          plus sling_pt*, sling_w*, sling_lineage on the object
##
## Two settings decide what this produces, and both are choices, not results:
##
##   OMEGA lets Slingshot return disjoint trees. Left at its default, the
##   minimum spanning tree attaches every cluster to one backbone, so
##   Osteoclast-like Mac is placed on a lineage despite PAGA connectivity
##   below 0.1 to every state at every neighbourhood size from 5 to 20.
##
##   START_CLUS is an assumption. Slingshot orders cells away from the cluster
##   you nominate; it cannot find the origin. Classical Mono is nominated from
##   Ccr2/Ly6c2 marker identity, and the methods section should say so.
##
## Resident Mac is excluded from fitting: 98.6% of those cells come from one
## library, and in the integrated embedding their distance to Remodeling Mac
## is shorter than Classical Mono to Transitional Mac, so no omega separates
## Resident Mac without also severing the Mono-Transitional edge that every
## within-library analysis recovered. Set DROP_STATES to character(0) to see
## this for yourself.
##
## Usage:
##   mono_mac <- readRDS("mono_mac_5state.rds")
##   source("slingshot_monomac.R")
##   saveRDS(mono_mac, "mono_mac_5state.rds")   # to keep the new columns
## =============================================================================

library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)


# ---- settings ---------------------------------------------------------------

OUT         <- "results_slingshot"
TAG         <- "noResident"

START_CLUS  <- "Classical Mono"          # ASSUMED, from marker identity
DROP_STATES <- c("Resident Mac")         # see the header
OMEGA       <- TRUE                      # allow disjoint trees
NPCS        <- 30
MIN_CELLS   <- 20                        # too few to place on a graph
W_CUT       <- 0.8                       # curve-weight cutoff for pseudotime
SEED        <- 42

STATE_LEVELS <- c("Classical Mono", "Transitional Mac", "Resident Mac",
                  "Remodeling Mac", "Osteoclast-like Mac")

COLS <- c("Classical Mono"      = "#D98C4A",
          "Transitional Mac"    = "#AFD8E2",
          "Resident Mac"        = "#7ED0F0",
          "Remodeling Mac"      = "#5A61A8",
          "Osteoclast-like Mac" = "#111455")

PT_SIZE   <- 0.9
FONT      <- "sans"                     # fontconfig default
CURVATURE <- -0.22                       # lineage arrow bend; 0 = straight

dir.create(OUT, showWarnings = FALSE)
set.seed(SEED)


# ---- 1. cells to fit --------------------------------------------------------

obj <- mono_mac[, !(mono_mac$myeloid_subtype_short %in% DROP_STATES)]

counts <- table(droplevels(factor(obj$myeloid_subtype_short)))
keep   <- names(counts)[counts >= MIN_CELLS]
sparse <- setdiff(names(counts), keep)

if (length(sparse)) {
  cat("below", MIN_CELLS, "cells, dropped:", paste(sparse, collapse = ", "),
      "\n  (too sparse to place, which is not the same as unconnected)\n")
}
obj <- obj[, obj$myeloid_subtype_short %in% keep]

state <- factor(as.character(obj$myeloid_subtype_short),
                levels = STATE_LEVELS[STATE_LEVELS %in% keep])

cat("\n", ncol(obj), " cells, ", nlevels(state), " states\n", sep = "")
print(table(state))

cat("\nlibrary composition per state (row fractions):\n")
print(round(prop.table(table(state, obj$sample), 1), 3))


# ---- 2. fit the trajectory --------------------------------------------------
# Fitted in harmony space. UMAP distorts distances non-monotonically, so a
# curve fitted on UMAP coordinates follows the picture, not the transcriptome.

sce <- SingleCellExperiment(
  assays = list(logcounts = GetAssayData(obj, assay = "RNA", layer = "data")))
reducedDims(sce) <- SimpleList(
  HARMONY = Embeddings(obj, "harmony")[, 1:NPCS],
  UMAP    = Embeddings(obj, "umap"))
colData(sce)$state <- as.character(state)

sce <- slingshot(sce,
                 clusterLabels = "state",
                 reducedDim    = "HARMONY",
                 start.clus    = START_CLUS,
                 omega         = OMEGA,
                 approx_points = 150)

sds      <- SlingshotDataSet(sce)
lineages <- slingLineages(sds)

cat("\n== lineages ==\n")
for (i in seq_along(lineages)) {
  cat(sprintf("  Lineage %d: %s\n", i, paste(lineages[[i]], collapse = " -> ")))
}

cat("\n== cluster MST edges (compare against the PAGA edge table) ==\n")
print(slingMST(sds))

orphans <- setdiff(levels(state), unique(unlist(lineages)))
if (length(orphans)) {
  cat("\non no lineage:", paste(orphans, collapse = ", "),
      "\n  With omega enabled this is a result, not a filtering step.\n")
}


# ---- 3. curve coordinates in UMAP space -------------------------------------
# embedCurves refits the same lineages and cell weights in the target
# embedding, so a drawn curve corresponds to the fitted one. Note that a curve
# fitted in 30 dimensions can self-intersect once projected to two, which is a
# property of the projection; the figures below therefore do not draw it.

curves <- slingCurves(embedCurves(sce, "UMAP"))
cat("\nembedded curves:", length(curves), "of", length(lineages), "lineages\n")

if (length(curves)) {
  curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
    pts <- curves[[i]]$s[curves[[i]]$ord, , drop = FALSE]
    data.frame(curve = i, order = seq_len(nrow(pts)),
               UMAP1 = pts[, 1], UMAP2 = pts[, 2])
  }))
  write.csv(curve_df, file.path(OUT, sprintf("curves_%s.csv", TAG)),
            row.names = FALSE)
}


# ---- 4. pseudotime, masked by curve weight ----------------------------------
# Slingshot gives every cell a projected pseudotime on every curve, including
# cells that never travel along it. Unmasked, a later FeaturePlot colours
# Osteoclast-like Mac somewhere in the middle of lineage 1. The masked columns
# are the ones to plot and report.

pt <- slingPseudotime(sce)
wt <- slingCurveWeights(sce)

put_on_object <- function(values, name) {
  v <- rep(NA_real_, ncol(mono_mac))
  names(v) <- colnames(mono_mac)
  v[colnames(sce)] <- values
  mono_mac[[name]] <- v
  mono_mac
}

for (i in seq_len(ncol(pt))) {
  masked <- pt[, i]
  masked[!(wt[, i] > W_CUT)] <- NA_real_
  mono_mac <- put_on_object(masked, sprintf("sling_pt%d", i))
  mono_mac <- put_on_object(wt[, i], sprintf("sling_w%d", i))
}

# Which lineage each cell belongs to: the curve it weighs most heavily on.
best <- apply(wt, 1, which.max)
best[apply(wt, 1, max) <= W_CUT] <- NA
mono_mac <- put_on_object(as.numeric(best), "sling_lineage")


# ---- 5. tables --------------------------------------------------------------

cat("\n== pseudotime quartiles by state, curve weight >", W_CUT, "==\n")
quart <- list()
for (i in seq_len(ncol(pt))) {
  on <- wt[, i] > W_CUT & is.finite(pt[, i])
  if (sum(on) < MIN_CELLS) {
    cat(sprintf("  curve %d: too few cells on it\n", i))
    next
  }
  q <- t(sapply(split(pt[on, i], droplevels(state[on])),
                quantile, c(0.25, 0.50, 0.75)))
  q <- data.frame(curve = i, state = rownames(q),
                  n = as.integer(table(droplevels(state[on]))[rownames(q)]),
                  round(q, 1), row.names = NULL, check.names = FALSE)
  cat(sprintf("  curve %d, %d cells\n", i, sum(on)))
  print(q, row.names = FALSE)
  quart[[length(quart) + 1]] <- q
}
if (length(quart)) {
  write.csv(do.call(rbind, quart),
            file.path(OUT, sprintf("pseudotime_quartiles_%s.csv", TAG)),
            row.names = FALSE)
}

cells <- data.frame(cell    = colnames(sce),
                    state   = as.character(state),
                    library = as.character(obj$sample),
                    UMAP1   = reducedDims(sce)$UMAP[, 1],
                    UMAP2   = reducedDims(sce)$UMAP[, 2],
                    pt, wt, check.names = FALSE)
write.csv(cells, file.path(OUT, sprintf("cells_%s.csv", TAG)), row.names = FALSE)


# ---- 6. figures -------------------------------------------------------------

# Plot from the object so the figures always match the stored columns.
plot_df <- data.frame(
  UMAP1 = Embeddings(mono_mac, "umap")[, 1],
  UMAP2 = Embeddings(mono_mac, "umap")[, 2],
  state = factor(as.character(mono_mac$myeloid_subtype_short),
                 levels = STATE_LEVELS),
  pt    = mono_mac$sling_pt1)

# Only lineage-1 states carry a pseudotime. Masking by weight alone leaves
# residual values on states that were never fitted, so restrict by state too.
lineage1 <- lineages[[1]]
plot_df$pt[!(plot_df$state %in% lineage1)] <- NA

cat("\nlineage 1:", paste(lineage1, collapse = " -> "), "\n")
cat(sum(is.finite(plot_df$pt)), "of", nrow(plot_df), "cells carry a pseudotime\n")

base_theme <- theme_void(base_family = FONT, base_size = 11) +
  theme(legend.position = "none",
        plot.margin = margin(8, 8, 8, 8))

# Lineage arrows run between cluster centroids rather than tracing the fitted
# curve. A curve fitted in 30 harmony dimensions self-intersects once
# projected to two UMAP dimensions, and that loop belongs to the projection,
# not to the trajectory. Centroid arrows carry the same lineage order without
# implying a path through UMAP space.
cent <- data.frame(
  state = levels(droplevels(plot_df$state)),
  UMAP1 = tapply(plot_df$UMAP1, droplevels(plot_df$state), median),
  UMAP2 = tapply(plot_df$UMAP2, droplevels(plot_df$state), median),
  row.names = NULL)

steps <- data.frame(from = head(lineage1, -1), to = tail(lineage1, -1))
steps$x    <- cent$UMAP1[match(steps$from, cent$state)]
steps$y    <- cent$UMAP2[match(steps$from, cent$state)]
steps$xend <- cent$UMAP1[match(steps$to,   cent$state)]
steps$yend <- cent$UMAP2[match(steps$to,   cent$state)]

# Pull both ends in so an arrow starts and stops outside the cluster cores.
TRIM <- 0.12
steps$x1 <- steps$x    + TRIM * (steps$xend - steps$x)
steps$y1 <- steps$y    + TRIM * (steps$yend - steps$y)
steps$x2 <- steps$xend - TRIM * (steps$xend - steps$x)
steps$y2 <- steps$yend - TRIM * (steps$yend - steps$y)

# White underlay first, so the arrow stays readable over dark clusters.
arrow_layers <- list(
  geom_curve(data = steps, aes(x = x1, y = y1, xend = x2, yend = y2),
             inherit.aes = FALSE, curvature = CURVATURE,
             colour = "white", linewidth = 1.7, lineend = "round",
             arrow = arrow(length = unit(3.3, "mm"), type = "open", angle = 25)),
  geom_curve(data = steps, aes(x = x1, y = y1, xend = x2, yend = y2),
             inherit.aes = FALSE, curvature = CURVATURE,
             colour = "#1A1A1A", linewidth = 0.85, lineend = "round",
             arrow = arrow(length = unit(3.0, "mm"), type = "open", angle = 25)))

p_states <- ggplot(plot_df, aes(UMAP1, UMAP2, colour = state)) +
  geom_point(size = PT_SIZE, stroke = 0) +
  scale_colour_manual(values = COLS, drop = FALSE) +
  arrow_layers +
  coord_equal() +
  base_theme

p_pseudotime <- ggplot(mapping = aes(UMAP1, UMAP2)) +
  geom_point(data = plot_df[is.na(plot_df$pt), ],
             colour = "#E4E4E4", size = PT_SIZE, stroke = 0) +
  geom_point(data = plot_df[!is.na(plot_df$pt), ], aes(colour = pt),
             size = PT_SIZE, stroke = 0) +
  scale_colour_viridis_c(name = "Pseudotime") +
  coord_equal() +
  base_theme +
  theme(legend.position = "right",
        legend.title = element_text(family = FONT, size = 9)) +
  guides(colour = guide_colourbar(barwidth = unit(3, "mm"),
                                  barheight = unit(28, "mm")))

ggsave(file.path(OUT, "fig_states.pdf"), p_states,
       width = 5.2, height = 5, device = cairo_pdf)
ggsave(file.path(OUT, "fig_pseudotime.pdf"), p_pseudotime,
       width = 5.6, height = 5, device = cairo_pdf)

if (requireNamespace("patchwork", quietly = TRUE)) {
  ggsave(file.path(OUT, "fig_combined.pdf"),
         patchwork::wrap_plots(p_states, p_pseudotime),
         width = 10.4, height = 5, device = cairo_pdf)
}


# ---- 7. provenance ----------------------------------------------------------

mono_mac@misc$slingshot <- list(
  tag           = TAG,
  lineages      = vapply(lineages, paste, character(1), collapse = " -> "),
  reduction     = sprintf("harmony (%d dims); curves embedded in umap", NPCS),
  start_cluster = paste(START_CLUS, "- ASSUMED from Ccr2/Ly6c2 identity"),
  dropped       = DROP_STATES,
  omega         = OMEGA,
  weight_cutoff = W_CUT,
  cells_fitted  = ncol(sce),
  date          = as.character(Sys.Date()))

mono_mac <<- mono_mac

cat("\nwrote", OUT, "/\n")
cat("added sling_pt*, sling_w*, sling_lineage to mono_mac\n")
cat("save with saveRDS(mono_mac, \"mono_mac_5state.rds\") to keep them\n")
