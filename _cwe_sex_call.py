#!/usr/bin/env python3
"""
Call the sex of each cell from RNA counts (AnnData), for pooled male/female samples.

Logic
  * Female signal  = XIST counts.
  * Male signal    = summed counts of Y-linked genes (non-PAR, so not shared with X).
  * Seeds: confidently female = XIST >= seed_xist; confidently male = Y counts >= seed_y.
           Cells meeting both are excluded from both seed sets (likely opposite-sex doublets).
  * Ambient levels are estimated from the seeds: the 99th percentile of XIST counts in the
    seed males and of Y counts in the seed females. These are what contamination from the
    other sex's RNA looks like in this dataset, so calls must exceed them.
  * female  : XIST > ambient_xist  and  Y <= ambient_y
    male    : Y > ambient_y        and  XIST <= ambient_xist
    ambiguous: both exceed  (doublets of opposite sex, or strong ambient RNA)
    unknown : neither exceeds (mostly dropout; expect many, especially for males)

Usage
    from sex_call_rna import call_sex, summarize_by_cluster, plot_diagnostic
    call_sex(a)                                   # adds a.obs["sex_call"] etc.
    tab = summarize_by_cluster(a, "atac_leiden")
    plot_diagnostic(a, "sex_diagnostic.png")
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "TXLNGY", "NLGN4Y", "TMSB4Y"]


def get_counts(adata, layer="counts"):
    """Return (matrix, gene_index) of raw integer counts, in order of preference."""
    if layer and layer in adata.layers:
        X, genes = adata.layers[layer], adata.var_names
    elif adata.raw is not None:
        X, genes = adata.raw.X, adata.raw.var_names
    else:
        X, genes = adata.X, adata.var_names
    data = X.data if sp.issparse(X) else np.asarray(X)
    sample = data[:10000]
    if len(sample) and not np.allclose(sample, np.round(sample)):
        raise ValueError("Counts look normalised/log-transformed. Point `layer=` at raw counts.")
    return X, pd.Index(genes)


def _sum_genes(X, genes, names):
    present = [g for g in names if g in genes]
    if not present:
        return np.zeros(X.shape[0]), np.zeros(X.shape[0]), present
    idx = [genes.get_loc(g) for g in present]
    sub = X[:, idx]
    sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    return sub.sum(axis=1), (sub > 0).sum(axis=1), present


def call_sex(adata, layer="counts", seed_xist=3, seed_y=2, y_genes=Y_GENES, q=0.99, verbose=True):
    X, genes = get_counts(adata, layer)
    if "XIST" not in genes:
        raise ValueError("XIST not found in var_names (are genes Ensembl IDs? map to symbols first).")
    xist = np.asarray(X[:, genes.get_loc("XIST")].todense() if sp.issparse(X)
                      else X[:, genes.get_loc("XIST")]).ravel()
    ysum, ygenes, present = _sum_genes(X, genes, y_genes)
    if verbose:
        print(f"Y genes found ({len(present)}): {present}")

    # Seeds use the positive signal only. Requiring "no Y" in females (or "no XIST" in males)
    # would force the ambient estimates to zero. Cells carrying both signals (opposite-sex
    # doublets) are excluded from both seed sets.
    both = (xist >= seed_xist) & (ysum >= seed_y)
    seed_f = (xist >= seed_xist) & ~both
    seed_m = (ysum >= seed_y) & ~both
    if seed_f.sum() == 0 or seed_m.sum() == 0:
        raise ValueError(f"Seeds empty (female={seed_f.sum()}, male={seed_m.sum()}); "
                         "lower seed_xist / seed_y or check the matrix.")
    amb_x = float(np.quantile(xist[seed_m], q))
    amb_y = float(np.quantile(ysum[seed_f], q))
    if verbose:
        print(f"seed cells: female={seed_f.sum()}, male={seed_m.sum()} | "
              f"ambient XIST (in males)={amb_x:.0f}, ambient Y (in females)={amb_y:.0f}")

    x_pos = xist > max(amb_x, 1)          # at least 2 XIST counts, and above male ambient level
    y_pos = ysum > amb_y
    call = np.select([x_pos & ~y_pos, y_pos & ~x_pos, x_pos & y_pos],
                     ["female", "male", "ambiguous"], default="unknown")
    adata.obs["xist_counts"] = xist
    adata.obs["y_counts"] = ysum
    adata.obs["y_genes_detected"] = ygenes
    adata.obs["sex_call"] = pd.Categorical(call, categories=["female", "male", "ambiguous", "unknown"])
    if verbose:
        print(adata.obs["sex_call"].value_counts())
    return adata


def summarize_by_cluster(adata, cluster_col):
    obs = adata.obs
    ct = pd.crosstab(obs[cluster_col], obs["sex_call"])
    sexed = ct[["female", "male"]].sum(axis=1)
    out = ct.copy()
    out["frac_female_of_sexed"] = (ct["female"] / sexed.replace(0, np.nan)).round(2)
    out["frac_unknown"] = (ct["unknown"] / ct.sum(axis=1)).round(2)
    tot = ct[["female", "male"]].sum()
    print(f"overall female fraction among sexed cells: {tot['female'] / tot.sum():.2f}")
    return out


def plot_diagnostic(adata, outfile="sex_diagnostic.png"):
    o = adata.obs
    colors = {"female": "tab:red", "male": "tab:blue", "ambiguous": "tab:purple", "unknown": "0.75"}
    rng = np.random.default_rng(0)
    jx = rng.uniform(-0.15, 0.15, len(o)); jy = rng.uniform(-0.15, 0.15, len(o))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for s, c in colors.items():
        m = (o["sex_call"] == s).to_numpy()
        ax[0].scatter(np.log1p(o["xist_counts"])[m] + jx[m], np.log1p(o["y_counts"])[m] + jy[m],
                      s=5, c=c, label=f"{s} ({m.sum()})", alpha=0.7)
    ax[0].set_xlabel("log1p(XIST counts)"); ax[0].set_ylabel("log1p(summed Y-gene counts)")
    ax[0].legend(markerscale=3, fontsize=8)
    ax[1].hist(np.log1p(o["xist_counts"]), bins=60, color="0.4")
    ax[1].set_xlabel("log1p(XIST counts)"); ax[1].set_ylabel("cells")
    ax[1].set_title("should be bimodal if both sexes are present")
    fig.tight_layout(); fig.savefig(outfile, dpi=150); plt.close(fig)
