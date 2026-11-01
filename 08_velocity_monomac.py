#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Integrated scVelo workflow for Li_AAA myeloid cells.

Stages
------
check : barcode matching between kb AnnData and Seurat metadata
fit   : matching, AnnData merge, dynamical RNA velocity, latent time; write h5ad
plot  : read saved h5ad and export publication-style PNG/PDF (no re-fitting)
all   : fit + plot (barcode checks are included in fit)

Usage
-----
python run_myeloid_RNA_velocity.py --stage check
python run_myeloid_RNA_velocity.py --stage fit --n-jobs 8
python run_myeloid_RNA_velocity.py --stage plot
python run_myeloid_RNA_velocity.py --stage all

Paths are relative to ~/Li_AAA by default; override with --root.
"""

import argparse
import gc
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  # Must be set before pyplot / scVelo plotting imports.
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import scvelo as scv


SAMPLES = {"AAA50": "50", "AAA100": "100", "AAA150": "150", "AAA200": "200"}
SUBTYPE_KEY = "myeloid_subtype_short"
CATEGORY_COLS = [
    "myeloid_subtype_short", "orig.ident", "group", "dose", "sample",
    "celltype_major", "myeloid_subtype", "myeloid_lineage",
]
SUBTYPE_COLORS = {
    "Classical Mono": "#D98C4A",
    "Transitional Mac": "#AFD8E2",
    "Resident Mac": "#7ED0F0",
    "Remodeling Mac": "#5A61A8",
    "Osteoclast-like Mac": "#111455",
    "Fibrocyte-like": "#8C6BB1",
    "Immature Neu": "#E7A9A1",
    "Inflammatory Neu": "#D4675F",
    "Mature/CXCR4hi Neu": "#A83E3E",
    "Migratory DC": "#65A98F",
    "cDC1": "#2F7F65",
    "Hb-high": "#BEBEBE",
    "Platelet-associated": "#8C8C8C",
}


def load_metadata(root):
    """Read Seurat metadata (must contain a 'cell' index)."""
    path = root / "05_velocity" / "meta.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing metadata: {path}")
    meta = pd.read_csv(path, index_col="cell")
    meta.index = meta.index.astype(str)
    if not meta.index.is_unique:
        raise ValueError("Duplicate cell barcodes in 05_velocity/meta.csv")
    if SUBTYPE_KEY not in meta.columns:
        raise KeyError(f"Required column not found in metadata: {SUBTYPE_KEY}")
    print(f"Seurat cells: {meta.shape[0]}", flush=True)
    return meta


def normalize_barcode(bc, sid):
    """Keep pre-existing -1 suffix if present; add sample prefix exactly once."""
    bc = str(bc)
    prefix = f"{sid}_"
    if bc.startswith(prefix):
        return bc
    if not bc.endswith("-1"):
        bc += "-1"
    return prefix + bc


def read_sample(root, folder, sid):
    p = root / f"velo_out_{folder}" / "counts_unfiltered"
    fn = p / "adata.h5ad"
    if not fn.exists():
        raise FileNotFoundError(f"Missing count AnnData: {fn}")
    a = ad.read_h5ad(fn)
    a.obs_names = pd.Index([normalize_barcode(bc, sid) for bc in a.obs_names])
    if not a.obs_names.is_unique:
        raise ValueError(f"Non-unique barcodes after prefixing in {folder}")
    if not a.var_names.is_unique:
        raise ValueError(f"Non-unique gene IDs in {folder}; cannot safely concat")
    return a, p


def check_barcodes(root, meta):
    """Stage 1, without constructing a merged dataset."""
    keep = set(meta.index)
    rows = []
    found = set()
    for folder, sid in SAMPLES.items():
        a, _ = read_sample(root, folder, sid)
        hit = a.obs_names.isin(keep)
        matches = a.obs_names[hit]
        rows.append({"sample": folder, "total_kb": a.n_obs, "recovered": len(matches),
                     "fraction": round(float(hit.mean()), 4)})
        print(f"{folder}: {len(matches)} recovered (kb barcodes: {a.n_obs})", flush=True)
        if len(matches) == 0:
            print("   kb example    :", list(a.obs_names[:2]))
            print("   Seurat example:", list(meta.index[meta.index.str.startswith(sid + "_")][:2]))
        found.update(matches)
        del a
    summary = pd.DataFrame(rows)
    outdir = root / "05_velocity"
    outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outdir / "barcode_recovery.csv", index=False)
    print(f"Total recovered: {len(found)} / {len(keep)} Seurat cells", flush=True)
    print(f"Saved: {outdir / 'barcode_recovery.csv'}", flush=True)
    return summary


def load_coordinate_table(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing coordinates: {path}")
    df = pd.read_csv(path, index_col="cell")
    df.index = df.index.astype(str)
    if not df.index.is_unique:
        raise ValueError(f"Duplicate barcodes in {path}")
    return df


def read_gene_symbols(path, n_genes):
    # The supplied kb file contains gene names (one per line).
    if not path.exists():
        raise FileNotFoundError(f"Missing gene names: {path}")
    symbols = pd.read_csv(path, header=None, sep="\t", dtype=str).iloc[:, 0]
    if len(symbols) != n_genes:
        raise ValueError(
            f"Gene name count mismatch in {path}: {len(symbols)} vs {n_genes} genes"
        )
    return symbols.fillna("").astype(str).tolist()


def build_anndata(root, meta):
    """Stage 2 input: align cells on barcode and genes on ORIGINAL gene IDs."""
    umap = load_coordinate_table(root / "05_velocity" / "umap.csv")
    harm = load_coordinate_table(root / "05_velocity" / "harmony.csv")
    if umap.shape[1] != 2:
        raise ValueError(f"Expected exactly 2 UMAP coordinates, got {umap.shape[1]}")

    keep = set(meta.index)
    adatas = []
    gene_to_symbol = {}
    recovered = set()
    recovery_rows = []

    for folder, sid in SAMPLES.items():
        a, p = read_sample(root, folder, sid)
        n_kb = a.n_obs
        matched = a.obs_names.isin(keep)
        a = a[matched].copy()
        print(f"{folder}: {a.n_obs}/{n_kb} matched cells; shape {a.shape}", flush=True)
        recovery_rows.append({"sample": folder, "total_kb": n_kb,
                              "recovered": a.n_obs,
                              "fraction": round(a.n_obs / n_kb, 4) if n_kb else 0})
        if a.n_obs == 0:
            raise ValueError(f"No matched cells for {folder}; run --stage check")
        if not {"mature", "ambiguous", "nascent"}.issubset(a.layers):
            missing = {"mature", "ambiguous", "nascent"} - set(a.layers)
            raise KeyError(f"Missing expected kb layers in {folder}: {sorted(missing)}")

        names = read_gene_symbols(p / "cells_x_genes.genes.names.txt", a.n_vars)
        for gene_id, symbol in zip(a.var_names.astype(str), names):
            symbol = symbol.strip() or str(gene_id)
            old = gene_to_symbol.get(gene_id)
            if old is not None and old != symbol:
                raise ValueError(
                    f"Conflicting gene symbols for ID {gene_id}: {old} vs {symbol}. "
                    "Check the annotation files before merging."
                )
            gene_to_symbol[gene_id] = symbol

        overlap = recovered.intersection(a.obs_names)
        if overlap:
            raise ValueError(f"Duplicate cells across samples, e.g. {next(iter(overlap))}")
        recovered.update(a.obs_names)
        adatas.append(a)

    pd.DataFrame(recovery_rows).to_csv(
        root / "05_velocity" / "barcode_recovery.csv", index=False
    )
    # Merge on stable gene IDs, not on column positions or on ambiguous symbols.
    # fill_value=0 also prevents NaNs from missing genes in the outer union.
    adata = ad.concat(adatas, join="outer", fill_value=0)
    gene_ids = adata.var_names.astype(str).tolist()
    adata.var["gene_id"] = gene_ids
    adata.var["symbol"] = [gene_to_symbol.get(gid, gid) for gid in gene_ids]
    adata.var_names = pd.Index(adata.var["symbol"].astype(str))
    adata.var_names_make_unique()
    del adatas
    gc.collect()

    missing_umap = adata.obs_names.difference(umap.index)
    missing_harm = adata.obs_names.difference(harm.index)
    if len(missing_umap) or len(missing_harm):
        raise ValueError(
            f"Missing coordinate rows: UMAP={len(missing_umap)}, Harmony={len(missing_harm)}"
        )

    # Exact barcode alignment; never rely on equal row order across CSV files.
    adata.obs = meta.loc[adata.obs_names].copy()
    adata.obsm["X_umap"] = umap.loc[adata.obs_names].to_numpy(dtype=np.float32)
    adata.obsm["X_harmony"] = harm.loc[adata.obs_names].to_numpy(dtype=np.float32)
    for col in CATEGORY_COLS:
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype("category")

    print(f"Merged: {adata.shape}; Harmony dims: {adata.obsm['X_harmony'].shape[1]}", flush=True)
    print(f"Seurat cells absent from kb input: {len(meta.index.difference(adata.obs_names))}", flush=True)
    return adata


def save_proportions(adata, outdir):
    if "orig.ident" not in adata.obs:
        print("WARNING: orig.ident absent; skipped proportions plot", flush=True)
        return
    outdir.mkdir(parents=True, exist_ok=True)
    scv.pl.proportions(adata, groupby="orig.ident", show=False)
    plt.gcf().savefig(outdir / "proportions.png", bbox_inches="tight", dpi=300)
    plt.close("all")


def velocity_path(root, ambiguous_mode):
    name = ("myeloid_velocity.h5ad" if ambiguous_mode == "spliced"
            else "myeloid_velocity_exclude_ambiguous.h5ad")
    return root / "05_velocity" / name


def fit_model(root, meta, n_jobs, ambiguous_mode):
    """Stage 2: dynamical scVelo analysis (the expensive stage)."""
    adata = build_anndata(root, meta)
    if ambiguous_mode == "spliced":
        # Keeps the original analysis assumption. Validate biologically with a
        # sensitivity run using --ambiguous-mode exclude.
        adata.layers["spliced"] = adata.layers["mature"] + adata.layers["ambiguous"]
    else:
        adata.layers["spliced"] = adata.layers["mature"].copy()
    adata.layers["unspliced"] = adata.layers["nascent"].copy()
    adata.uns["ambiguous_mode"] = ambiguous_mode

    save_proportions(adata, root / "05_velocity" / "figures")

    print("Filter/normalize...", flush=True)
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    print("Moments on Harmony coordinates...", flush=True)
    scv.pp.moments(adata, n_pcs=None, n_neighbors=30, use_rep="X_harmony")
    print("Recover dynamics...", flush=True)
    scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata, n_jobs=n_jobs)
    scv.tl.velocity_confidence(adata)
    scv.tl.latent_time(adata)

    output = velocity_path(root, ambiguous_mode)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(f"Saved model: {output}", flush=True)
    return adata


def configure_plotting():
    have = {f.name for f in fm.fontManager.ttflist}
    preferred = ["Arial", "Liberation Sans", "Nimbus Sans", "FreeSans", "DejaVu Sans"]
    chosen = next((name for name in preferred if name in have), "DejaVu Sans")
    print(f"Font in use: {chosen}", flush=True)
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = [chosen] + preferred
    matplotlib.rcParams["axes.unicode_minus"] = False
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["savefig.dpi"] = 600
    matplotlib.rcParams["figure.dpi"] = 150


def set_subtype_palette(adata):
    if SUBTYPE_KEY not in adata.obs:
        raise KeyError(f"Missing {SUBTYPE_KEY} in saved AnnData")
    observed = set(adata.obs[SUBTYPE_KEY].dropna().astype(str))
    order = [name for name in SUBTYPE_COLORS if name in observed]
    unknown = sorted(observed.difference(SUBTYPE_COLORS))
    if unknown:
        print(f"WARNING: no preset color for {unknown}; using gray", flush=True)
    order.extend(unknown)
    adata.obs[SUBTYPE_KEY] = pd.Categorical(
        adata.obs[SUBTYPE_KEY], categories=order, ordered=True
    )
    adata.uns[f"{SUBTYPE_KEY}_colors"] = [SUBTYPE_COLORS.get(c, "#BBBBBB") for c in order]


def save_figure(fig, outdir, basename):
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / (basename + ".png")
    pdf = outdir / (basename + ".pdf")
    fig.savefig(png, bbox_inches="tight", dpi=600)
    try:
        fig.savefig(pdf, bbox_inches="tight")
        print(f"Saved {basename}: PNG + PDF", flush=True)
    except (ValueError, RuntimeError) as exc:
        print(f"PDF vector export failed for {basename}: {exc}", flush=True)
        # Rasterize collection artists only; keep labels as vector text.
        try:
            for axis in fig.axes:
                for artist in axis.collections:
                    artist.set_rasterized(True)
            fig.savefig(pdf, bbox_inches="tight", dpi=600)
            print(f"Saved {basename}: PNG + PDF (collections rasterized)", flush=True)
        except Exception as err:
            print(f"WARNING: PDF not available for {basename}: {err}", flush=True)
    finally:
        plt.close(fig)


def make_plots(root, adata=None, ambiguous_mode="spliced"):
    """Stage 3: publication plots without re-fitting velocity model."""
    if adata is None:
        path = velocity_path(root, ambiguous_mode)
        if not path.exists():
            raise FileNotFoundError(f"Missing fitted velocity result: {path}. Run --stage fit")
        adata = scv.read(str(path))

    configure_plotting()
    set_subtype_palette(adata)
    folder = ("figures_final" if ambiguous_mode == "spliced"
              else "figures_final_exclude_ambiguous")
    outdir = root / "05_velocity" / folder
    if "orig.ident" not in adata.obs:
        raise KeyError("orig.ident is missing; needed for stream_group")
    for field in ("velocity_confidence", "latent_time"):
        if field not in adata.obs:
            raise KeyError(f"Missing analysis output: {field}. Run --stage fit")

    ax = scv.pl.velocity_embedding_stream(
        adata, basis="umap", color=SUBTYPE_KEY, legend_loc="right",
        title="RNA velocity — myeloid subtypes", show=False
    )
    save_figure(ax.get_figure(), outdir, "stream_subtype")

    ax = scv.pl.velocity_embedding_stream(
        adata, basis="umap", color="orig.ident", legend_loc="right", show=False
    )
    save_figure(ax.get_figure(), outdir, "stream_group")

    ax = scv.pl.scatter(adata, color="velocity_confidence", cmap="coolwarm", show=False)
    save_figure(ax.get_figure(), outdir, "confidence")

    ax = scv.pl.scatter(adata, color="latent_time", color_map="gnuplot", show=False)
    save_figure(ax.get_figure(), outdir, "latent_time")

    print("\n=== velocity_confidence by subtype ===", flush=True)
    summary = adata.obs.groupby(SUBTYPE_KEY, observed=True)["velocity_confidence"].agg(
        ["mean", "median", "count"]
    )
    print(summary.to_string(), flush=True)
    summary.to_csv(outdir / "confidence_by_subtype.csv")
    if "fit_alpha" in adata.var:
        print(f"Fitted genes: {adata.var['fit_alpha'].notna().sum()} / {adata.n_vars}", flush=True)
    print(f"Final figures: {outdir}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("~/Li_AAA"),
                        help="Li_AAA project root (default: ~/Li_AAA)")
    parser.add_argument("--stage", choices=["check", "fit", "plot", "all"],
                        default="all")
    parser.add_argument("--n-jobs", type=int, default=8)
    parser.add_argument("--ambiguous-mode", choices=["spliced", "exclude"],
                        default="spliced", help="Original model: mature+ambiguous as spliced")
    args = parser.parse_args()
    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Project folder does not exist: {root}")
    print(f"Project: {root}; stage: {args.stage}", flush=True)

    if args.stage == "plot":
        make_plots(root, ambiguous_mode=args.ambiguous_mode)
    else:
        meta = load_metadata(root)
        if args.stage == "check":
            check_barcodes(root, meta)
        else:
            # Fit includes barcode matching, so --stage all does not reread inputs.
            adata = fit_model(root, meta, args.n_jobs, args.ambiguous_mode)
            if args.stage == "all":
                make_plots(root, adata=adata, ambiguous_mode=args.ambiguous_mode)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()

