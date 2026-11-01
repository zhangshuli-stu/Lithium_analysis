#!/usr/bin/env python
"""
=============================================================================
PAGA connectivity of the monocyte-macrophage compartment

Input :  paga_export/           written by export_monomac.R
Output:  results_paga/
           sample<L>/edges.csv         connectivity, all state pairs
           sample<L>/sensitivity.csv   connectivity vs n_neighbors and n_pcs
           sample<L>/knn_<state>.csv   sweep for states below SMALL_STATE
           sample<L>/paga.pdf          graph figure
           sample<L>/<L>.h5ad
           merged/...                  same, on the harmony embedding

Which run is the evidence
  In this design sample, dose and library are one variable, so a merged run on
  harmony cannot separate biology from integration. A state confined to one
  library gets pulled towards whichever cluster it most resembles elsewhere,
  and the resulting edge is unfalsifiable. The within-library runs carry no
  batch structure at all and are the evidence; the merged run is an overview.
  An edge that appears only in the merged graph is an integration result.

Comparing across libraries
  Connectivity is observed edges over edges expected at random, and the
  expectation depends on how many clusters are in the graph and how large each
  is. The four libraries differ enormously in composition, so their raw values
  are not on a common scale. Compare the ORDER of edges within a library
  across libraries, not the numbers.

Usage:
  python paga_monomac.py --mode within           # per library, the evidence
  python paga_monomac.py --mode merged           # harmony, the overview
  python paga_monomac.py --mode both
=============================================================================
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import FancyArrowPatch


# ---- settings ---------------------------------------------------------------

INDIR = "paga_export"
OUT = "results_paga"

STATE_LEVELS = ["Classical Mono", "Transitional Mac", "Resident Mac",
                "Remodeling Mac", "Osteoclast-like Mac"]

COLS = {"Classical Mono":      "#D98C4A",
        "Transitional Mac":    "#AFD8E2",
        "Resident Mac":        "#7ED0F0",
        "Remodeling Mac":      "#5A61A8",
        "Osteoclast-like Mac": "#111455"}

MIN_CELLS = 20          # too few to place on a graph
SMALL_STATE = 100       # below this, run the neighbourhood sweep
N_NEIGHBORS = 15
N_PCS = 30              # harmony has 50; the later ones are mostly noise
SEED = 0

# Thresholds for the wording only. 0.296 and 0.30 are the same number; what
# carries meaning is the ordering of edges and the gap between groups of them.
STRONG, WEAK = 0.30, 0.10

# The topology under test, written as edges that must hold and edges that must
# not. Stating the falsifiers up front is what stops the result being read
# back into whichever shape it happens to take.
HYPOTHESIS = [("Classical Mono", "Transitional Mac"),
              ("Transitional Mac", "Resident Mac"),
              ("Transitional Mac", "Remodeling Mac"),
              ("Remodeling Mac", "Osteoclast-like Mac")]
FALSIFIERS = [("Classical Mono", "Resident Mac"),
              ("Classical Mono", "Remodeling Mac"),
              ("Classical Mono", "Osteoclast-like Mac"),
              ("Transitional Mac", "Osteoclast-like Mac"),
              ("Resident Mac", "Osteoclast-like Mac"),
              ("Resident Mac", "Remodeling Mac")]

THRESHOLDS = [0.05, 0.10, 0.20, 0.30]   # one figure per threshold
PLOT_THRESHOLD = 0.10
LABEL_SIZE = 13
FONT = "Arial"                          # falls back to Liberation Sans


# ---- helpers ----------------------------------------------------------------

def pick_font(preferred=FONT):
    avail = {f.name for f in fm.fontManager.ttflist}
    for c in (preferred, "Liberation Sans", "Nimbus Sans", "Helvetica",
              "DejaVu Sans"):
        if c in avail:
            return c
    return "sans-serif"


def load_h5ad(indir):
    """Assemble the flat files from export_monomac.R into an AnnData."""
    cfg = {}
    p = os.path.join(indir, "export_config.txt")
    if os.path.exists(p):
        for line in open(p):
            if "=" in line:
                k, v = line.strip().split("=", 1)
                cfg[k] = v

    X = scipy.io.mmread(os.path.join(indir, "X_data.mtx")).tocsr()
    genes = pd.read_csv(os.path.join(indir, "genes.csv"))["gene"].astype(str)
    cells = pd.read_csv(os.path.join(indir, "cells.csv"))["cell"].astype(str)
    obs = pd.read_csv(os.path.join(indir, "metadata.csv"), index_col=0)
    emb = pd.read_csv(os.path.join(indir, "embedding.csv"), index_col=0)

    adata = sc.AnnData(X=X.astype(np.float32),
                       obs=pd.DataFrame(index=cells),
                       var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()

    # Align on barcode; row order is not to be trusted across five files.
    assert set(cells) == set(obs.index), "barcode mismatch: metadata"
    assert set(cells) == set(emb.index), "barcode mismatch: embedding"
    adata.obs = obs.loc[adata.obs_names].copy()

    red = cfg.get("reduction", "harmony")
    adata.obsm[f"X_{red}"] = emb.loc[adata.obs_names].to_numpy(np.float32)

    up = os.path.join(indir, "umap.csv")
    if os.path.exists(up):
        umap = pd.read_csv(up, index_col=0)
        adata.obsm["X_umap"] = umap.loc[adata.obs_names].to_numpy(np.float32)

    # .X is already lognormalised; do not normalise again.
    adata.uns["config"] = cfg
    return adata, cfg


def graph_within(adata, key, n_neighbors, n_pcs, seed):
    """HVG, scale and PCA inside this library only.

    Features selected on the full object are ranked partly by between-library
    variance, which is exactly what a within-library run is meant to exclude.
    """
    sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars - 1))
    hv = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=min(50, hv.n_obs - 1, hv.n_vars - 1),
              svd_solver="arpack", random_state=seed)
    adata.obsm["X_pca"] = hv.obsm["X_pca"]
    return graph_from_rep(adata, key, "X_pca", n_neighbors, n_pcs, seed)


def graph_from_rep(adata, key, rep, n_neighbors, n_pcs, seed):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=rep,
                    n_pcs=min(n_pcs, adata.obsm[rep].shape[1]),
                    random_state=seed)
    sc.tl.paga(adata, groups=key)
    return adata


def edge_table(adata, key):
    cats = list(adata.obs[key].cat.categories)
    C = adata.uns["paga"]["connectivities"].toarray()
    rows = [{"state_a": cats[i], "state_b": cats[j], "connectivity": C[i, j]}
            for i, j in itertools.combinations(range(len(cats)), 2)]
    return pd.DataFrame(rows).sort_values("connectivity", ascending=False) \
                             .reset_index(drop=True)


def lookup(df, a, b):
    m = ((df.state_a == a) & (df.state_b == b)) | \
        ((df.state_a == b) & (df.state_b == a))
    return float(df.loc[m, "connectivity"].iloc[0]) if m.any() else np.nan


def verdict(v):
    if not np.isfinite(v):
        return "n/a"
    return "PRESENT" if v >= STRONG else "absent" if v <= WEAK else "ambiguous"


def knn_composition(adata, key, state):
    """Where the neighbours of one state come from.

    Connectivity near zero for a small cluster has two readings: the state is
    discrete, or its kNN budget is spent inside itself so no cross-cluster
    edge can form. This tells them apart.
    """
    from collections import Counter
    D = adata.obsp["distances"].tolil()
    labs = adata.obs[key].astype(str).values
    c = Counter()
    for i in np.flatnonzero(labs == state):
        for j in D.rows[i]:
            c[labs[j]] += 1
    total = max(sum(c.values()), 1)
    return pd.Series({k: v / total for k, v in c.items()}).sort_values(
        ascending=False)


# ---- figure -----------------------------------------------------------------

def layout(conn, threshold, seed=0):
    import igraph as ig
    n = conn.shape[0]
    g, w = ig.Graph(n=n), []
    for i in range(n):
        for j in range(i + 1, n):
            if conn[i, j] > threshold:
                g.add_edge(i, j)
                w.append(conn[i, j])
    if g.ecount() == 0:
        a = np.linspace(0, 2 * np.pi, n, endpoint=False)
        return np.column_stack([np.cos(a), np.sin(a)])
    pos = np.array(g.layout_fruchterman_reingold(weights=w).coords)
    span = pos.max(0) - pos.min(0)
    span[span == 0] = 1
    return (pos - pos.min(0)) / span


def plot_paga(adata, key, threshold, out, pos=None, edge_labels=False,
              label_size=LABEL_SIZE, label_pad=0.05, seed=SEED):
    cats = list(adata.obs[key].cat.categories)
    conn = adata.uns["paga"]["connectivities"].toarray()
    counts = adata.obs[key].value_counts()[cats].values

    plt.rcParams.update({"font.family": pick_font(), "pdf.fonttype": 42,
                         "ps.fonttype": 42})
    if pos is None:
        pos = layout(conn, threshold, seed)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))

    cmax = max(conn.max(), 1e-9)
    for i, j in itertools.combinations(range(len(cats)), 2):
        c = conn[i, j]
        if c <= threshold:
            continue
        ax.plot(*zip(pos[i], pos[j]), color="#333333",
                lw=0.4 + 1.9 * (c / cmax), solid_capstyle="round",
                zorder=1, alpha=0.85)
        if edge_labels:
            m = (pos[i] + pos[j]) / 2
            ax.text(m[0], m[1], f"{c:.2f}", fontsize=label_size * 0.7,
                    color="#555555", ha="center", va="center", zorder=2,
                    bbox=dict(fc="white", ec="none", pad=0.8, alpha=0.85))

    sizes = 50 + 260 * (counts / counts.max())
    ax.scatter(pos[:, 0], pos[:, 1], s=sizes,
               c=[COLS.get(c, "#CCCCCC") for c in cats],
               edgecolors="white", linewidths=0.9, zorder=3)

    # Labels are pushed out from the graph centre, measured from the node edge
    # so large and small nodes end up with the same visual gap.
    centre = pos.mean(axis=0)
    spread = np.linalg.norm(pos - centre, axis=1).max() or 1.0
    span = max(pos[:, 0].ptp(), pos[:, 1].ptp()) or 1.0
    radii = np.sqrt(sizes / np.pi) / 72.0 / 5.2 * span * 1.6
    for i, name in enumerate(cats):
        v = pos[i] - centre
        d = np.linalg.norm(v)
        v = v / d if d > 0.15 * spread else np.array([0.0, -1.0])
        p = pos[i] + v * (radii[i] + label_pad)
        ax.text(p[0], p[1], name, fontsize=label_size,
                ha="left" if v[0] > 0.25 else "right" if v[0] < -0.25 else "center",
                va="bottom" if v[1] > 0.25 else "top" if v[1] < -0.25 else "center",
                zorder=4)

    ax.set_xlim(pos[:, 0].min() - 0.45, pos[:, 0].max() + 0.45)
    ax.set_ylim(pos[:, 1].min() - 0.35, pos[:, 1].max() + 0.35)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return pos


# ---- one run ----------------------------------------------------------------

def run(adata, key, tag, outdir, rep=None, n_neighbors=N_NEIGHBORS,
        n_pcs=N_PCS, seed=SEED):
    os.makedirs(outdir, exist_ok=True)
    counts = adata.obs[key].value_counts()

    keep = counts[counts >= MIN_CELLS].index.tolist()
    dropped = [s for s in counts.index if s not in keep]
    if dropped:
        print(f"  below {MIN_CELLS} cells, dropped: {', '.join(dropped)}")
        print("  (too sparse to place, which is not the same as unconnected)")
    adata = adata[adata.obs[key].isin(keep)].copy()
    adata.obs[key] = adata.obs[key].cat.remove_unused_categories()

    # A state near the floor cannot support the default neighbourhood.
    nmin = int(adata.obs[key].value_counts().min())
    nn = min(n_neighbors, max(5, nmin // 3))
    if nn != n_neighbors:
        print(f"  n_neighbors {n_neighbors} -> {nn} "
              f"(smallest state has {nmin} cells)")

    adata = (graph_from_rep(adata, key, rep, nn, n_pcs, seed) if rep
             else graph_within(adata, key, nn, n_pcs, seed))

    edges = edge_table(adata, key)
    edges["n_a"] = edges.state_a.map(counts)
    edges["n_b"] = edges.state_b.map(counts)
    edges["call"] = edges.connectivity.map(verdict)
    edges.to_csv(f"{outdir}/edges.csv", index=False)
    print(f"\n== {tag}: connectivity ==")
    print(edges.to_string(index=False))

    print("\n  hypothesised edges")
    for a, b in HYPOTHESIS:
        if a in keep and b in keep:
            v = lookup(edges, a, b)
            print(f"    {a:22s} -- {b:22s} {v:.3f}  {verdict(v)}")
        else:
            print(f"    {a:22s} -- {b:22s}   not testable here")

    print("\n  falsifying edges")
    for a, b in FALSIFIERS:
        if a in keep and b in keep:
            v = lookup(edges, a, b)
            flag = "  <- contradicts the topology" if v >= STRONG else ""
            print(f"    {a:22s} -- {b:22s} {v:.3f}  {verdict(v)}{flag}")

    # Small states: is the isolation real, or is the neighbourhood too big?
    for st in [s for s in keep if counts[s] < SMALL_STATE]:
        print(f"\n  kNN composition of {st} (n={counts[st]}, n_neighbors={nn})")
        comp = knn_composition(adata, key, st)
        print("   " + comp.round(3).to_string().replace("\n", "\n   "))

        rows = []
        for k in (5, 8, 10, 15, 20):
            tmp = (graph_from_rep(adata.copy(), key, rep, k, n_pcs, seed) if rep
                   else graph_within(adata.copy(), key, k, n_pcs, seed))
            d = edge_table(tmp, key)
            r = {"n_neighbors": k,
                 "self_frac": knn_composition(tmp, key, st).get(st, np.nan)}
            for other in keep:
                if other != st:
                    r[other] = lookup(d, st, other)
            rows.append(r)
        sweep = pd.DataFrame(rows)
        sweep.to_csv(f"{outdir}/knn_{st.replace(' ', '_')}.csv", index=False)
        print(f"\n  {st}: connectivity vs n_neighbors")
        print("   " + sweep.round(3).to_string(index=False).replace("\n", "\n   "))

        best = sweep.drop(columns=["n_neighbors", "self_frac"]).max().max()
        if best < WEAK:
            print(f"\n    Highest connectivity at any neighbourhood size:"
                  f" {best:.3f}. This state does not attach to the graph.")

    # Sensitivity. An edge that appears at one setting only is not an edge.
    rows = []
    for k in (10, 15, 30, 50):
        tmp = (graph_from_rep(adata.copy(), key, rep, k, n_pcs, seed) if rep
               else graph_within(adata.copy(), key, k, n_pcs, seed))
        d = edge_table(tmp, key)
        for a, b in itertools.combinations(keep, 2):
            rows.append({"parameter": "n_neighbors", "value": k,
                         "edge": f"{a} -- {b}", "connectivity": lookup(d, a, b)})
    for p in (20, 30, 50):
        tmp = (graph_from_rep(adata.copy(), key, rep, nn, p, seed) if rep
               else graph_within(adata.copy(), key, nn, p, seed))
        d = edge_table(tmp, key)
        for a, b in itertools.combinations(keep, 2):
            rows.append({"parameter": "n_pcs", "value": p,
                         "edge": f"{a} -- {b}", "connectivity": lookup(d, a, b)})
    sens = pd.DataFrame(rows)
    sens.to_csv(f"{outdir}/sensitivity.csv", index=False)
    print(f"\n  sensitivity (written to {outdir}/sensitivity.csv)")
    print("   " + sens.pivot_table(index="edge", columns=["parameter", "value"],
                                   values="connectivity").round(3)
          .to_string().replace("\n", "\n   "))

    conn = adata.uns["paga"]["connectivities"].toarray()
    for thr in THRESHOLDS:
        if (conn > thr).sum() == 0:
            print(f"  (no edges above {thr}; figure skipped)")
            continue
        plot_paga(adata, key, thr, f"{outdir}/paga_thr{thr}.pdf",
                  edge_labels=(thr == PLOT_THRESHOLD))
    adata.write(f"{outdir}/{tag}.h5ad")
    print(f"\n  wrote {outdir}/")


# ---- main -------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default=INDIR)
    ap.add_argument("--outdir", default=OUT)
    ap.add_argument("--mode", default="both",
                    choices=["within", "merged", "both"])
    ap.add_argument("--n-neighbors", type=int, default=N_NEIGHBORS)
    ap.add_argument("--n-pcs", type=int, default=N_PCS)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sc.settings.verbosity = 1

    adata, cfg = load_h5ad(args.indir)
    key = cfg.get("label_col", "myeloid_subtype_short")
    skey = cfg.get("sample_col", "sample")
    rep = f"X_{cfg.get('reduction', 'harmony')}"

    adata.obs[key] = pd.Categorical(
        adata.obs[key], categories=[s for s in STATE_LEVELS
                                    if s in set(adata.obs[key])])
    print(f"{adata.n_obs} cells, {adata.obs[key].nunique()} states")
    print(adata.obs[key].value_counts().to_string())

    comp = pd.crosstab(adata.obs[key], adata.obs[skey])
    frac = comp.div(comp.sum(axis=1), axis=0)
    print("\nlibrary composition per state (row fractions)")
    print(frac.round(3).to_string())
    comp.to_csv(f"{args.outdir}/state_by_library.csv")
    for st, f in frac.max(axis=1)[frac.max(axis=1) > 0.8].items():
        print(f"  !! {st}: {f:.1%} from library {frac.loc[st].idxmax()}; "
              "every edge touching it is partly about that library")

    if args.mode in ("within", "both"):
        for lib in sorted(adata.obs[skey].astype(str).unique()):
            sub = adata[adata.obs[skey].astype(str) == lib].copy()
            sub.obs[key] = sub.obs[key].cat.remove_unused_categories()
            print(f"\n{'=' * 60}\nlibrary {lib}: {sub.n_obs} cells")
            print(sub.obs[key].value_counts().to_string())
            run(sub, key, lib, f"{args.outdir}/sample{lib}",
                rep=None, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs,
                seed=args.seed)

    if args.mode in ("merged", "both"):
        print(f"\n{'=' * 60}\nmerged on {rep}: {adata.n_obs} cells")
        print("Overview only. Dose and library are one variable here, so an")
        print("edge seen only in this run is an integration result.")
        run(adata.copy(), key, "merged", f"{args.outdir}/merged",
            rep=rep, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs,
            seed=args.seed)


if __name__ == "__main__":
    main()
