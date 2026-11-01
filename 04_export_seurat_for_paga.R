## =============================================================================
## Export mono_mac for PAGA
##
## Input :  mono_mac_5state.rds
## Output:  paga_export/
##            embedding.csv   harmony coordinates, what PAGA actually runs on
##            metadata.csv    state, library, QC columns
##            umap.csv        the existing UMAP, so figures line up
##            X_data.mtx      lognormalised expression, cells x genes
##            genes.csv, cells.csv, export_config.txt
##
## Usage:
##   mono_mac <- readRDS("mono_mac_5state.rds")
##   source("export_monomac.R")
## =============================================================================

library(Seurat)
library(Matrix)


# ---- settings ---------------------------------------------------------------

OUT        <- "paga_export"
LABEL_COL  <- "myeloid_subtype_short"
SAMPLE_COL <- "sample"
REDUCTION  <- "harmony"
UMAP_RED   <- "umap"

# Kept beyond the variable features so marker dynamics can be plotted along
# the trajectory without re-exporting.
EXTRA_GENES <- c("Ccr2", "Ly6c2", "Plac8", "Sell",
                 "Cd163", "Selenop", "Mrc1", "Gas6", "Folr2", "Lyve1",
                 "Gdf15", "Ccl8",
                 "Spp1", "Mmp12", "Mmp9", "Trem2", "Gpnmb", "Lgals3",
                 "Ctsk", "Atp6v0d2", "Acp5", "Nfatc1", "Dcstamp", "Oscar",
                 "Mki67", "Adgre1", "Itgam")

dir.create(OUT, showWarnings = FALSE)


# ---- 1. checks --------------------------------------------------------------

stopifnot(REDUCTION %in% Reductions(mono_mac),
          UMAP_RED  %in% Reductions(mono_mac),
          LABEL_COL %in% colnames(mono_mac@meta.data))

cat(ncol(mono_mac), "cells\n")
print(table(mono_mac@meta.data[[LABEL_COL]]))

cat("\nstate x library:\n")
print(table(mono_mac@meta.data[[LABEL_COL]], mono_mac@meta.data[[SAMPLE_COL]]))

cat("\nrow fractions:\n")
print(round(prop.table(table(mono_mac@meta.data[[LABEL_COL]],
                             mono_mac@meta.data[[SAMPLE_COL]]), 1), 3))
cat("A state confined to one library carries that library into every edge it\n",
    "appears in, which is why the Python script treats the merged run as an\n",
    "overview and the within-library runs as the evidence.\n", sep = "")


# ---- 2. export --------------------------------------------------------------

obj <- mono_mac
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

emb <- Embeddings(obj, REDUCTION)
write.csv(emb, file.path(OUT, "embedding.csv"), quote = FALSE)
cat(sprintf("\nembedding: %d x %d\n", nrow(emb), ncol(emb)))

write.csv(Embeddings(obj, UMAP_RED), file.path(OUT, "umap.csv"), quote = FALSE)

md <- obj@meta.data
md[[LABEL_COL]] <- as.character(md[[LABEL_COL]])
write.csv(md, file.path(OUT, "metadata.csv"), quote = TRUE)

genes <- union(VariableFeatures(obj), EXTRA_GENES)
genes <- intersect(genes, rownames(obj))
cat("genes exported:", length(genes), "\n")

# .X is Seurat's data slot: lognormalised, not scaled. The Python side must
# not normalise it again.
X <- GetAssayData(obj, assay = "RNA", layer = "data")[genes, ]
writeMM(as(t(X), "CsparseMatrix"), file.path(OUT, "X_data.mtx"))
write.csv(data.frame(gene = rownames(X)), file.path(OUT, "genes.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(data.frame(cell = colnames(X)), file.path(OUT, "cells.csv"),
          row.names = FALSE, quote = FALSE)

writeLines(c(paste0("label_col=", LABEL_COL),
             paste0("sample_col=", SAMPLE_COL),
             paste0("reduction=", REDUCTION),
             paste0("n_dims=", ncol(emb)),
             paste0("n_cells=", ncol(obj)),
             paste0("expression=seurat_data_lognorm"),
             paste0("date=", Sys.Date())),
           file.path(OUT, "export_config.txt"))

cat("\nwrote", OUT, "/\n")
