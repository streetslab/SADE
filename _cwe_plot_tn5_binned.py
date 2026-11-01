#!/usr/bin/env python3
"""
Per-cluster Tn5 insertion profiles for ONE chromosome whose insertions are already
binned into fixed windows (default 3 kb).

Assumed inputs
    insertions : dict  barcode -> counts per window for the loaded chromosome.
                 Accepted per-barcode values: 1-D numpy array, scipy sparse vector/matrix,
                 or dict {window_index: count}.  (Edit to_vector() if yours differs.)
    clusters   : pandas Series  barcode -> cluster label

Method
  * Each cell's window counts are scaled to sum to 1 (equal weight per cell), summed per
    cluster and renormalised -> cluster profile p_c.
  * Background = all cells pooled, equal weight per cell.
  * log2 enrichment = log2((p_c + pseudo) / (p_bg + pseudo)).
  * A region is called only if log2 enrichment >= min_log2 AND at least min_cell_frac of the
    cluster's cells have an insertion in it, so a few deep cells cannot create a "peak".
  * `rebin` merges k adjacent windows (k=10 -> 30 kb) for plotting/calling; detection is
    computed after merging, per cell.
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def to_vector(entry, n_bins):
    """Return a dense float vector of length n_bins. EDIT if your format differs."""
    if entry is None:                                   # stored as None = no insertions
        return np.zeros(n_bins)
    if isinstance(entry, dict):
        v = np.zeros(n_bins)
        idx = np.fromiter(entry.keys(), dtype=int, count=len(entry))
        v[idx] = np.fromiter(entry.values(), dtype=float, count=len(entry))
        return v
    if hasattr(entry, "rows") and entry.shape[0] == 1:  # lil_matrix (1 x n): fast path
        v = np.zeros(n_bins)
        cols = np.asarray(entry.rows[0], dtype=int)
        vals = np.array([0.0 if x is None else x for x in entry.data[0]], dtype=float)
        keep = cols < n_bins
        v[cols[keep]] = vals[keep]
        return v
    if hasattr(entry, "toarray"):                       # other scipy sparse
        entry = entry.toarray()
    v = np.asarray(entry, dtype=float).ravel()
    if len(v) < n_bins:
        v = np.pad(v, (0, n_bins - len(v)))
    return v[:n_bins]


def build_profile(insertions, clusters, n_bins=None, bin_size=3000, rebin=1):
    # use .get() so a defaultdict does not create new keys; barcodes with no entry or None
    # (no insertions on this chromosome) contribute nothing and are skipped below
    if n_bins is None:                                   # infer from the first non-empty entry
        first = next((insertions.get(bc) for bc in clusters.index
                      if insertions.get(bc) is not None), None)
        if first is None:
            raise ValueError("No cell has insertions on this chromosome; pass n_bins explicitly.")
        if isinstance(first, dict):
            n_bins = max(first) + 1          # may undercount; pass n_bins= if the last windows are empty
        elif hasattr(first, "shape"):
            n_bins = max(first.shape)        # dense 1-D or sparse (1 x n) / (n x 1)
        else:
            n_bins = len(first)
    n_out = int(np.ceil(n_bins / rebin))
    labels = sorted(clusters.unique(), key=str)
    prof = {k: np.zeros(n_out) for k in labels}
    det = {k: np.zeros(n_out) for k in labels}
    n_cells = {k: 0 for k in labels}
    for bc, k in clusters.items():
        entry = insertions.get(bc)
        if entry is None:
            continue
        v = to_vector(entry, n_bins)
        if rebin > 1:
            v = np.pad(v, (0, n_out * rebin - len(v))).reshape(n_out, rebin).sum(axis=1)
        s = v.sum()
        if s == 0:
            continue
        prof[k] += v / s
        det[k] += v > 0
        n_cells[k] += 1
    P = pd.DataFrame(prof).T
    D = pd.DataFrame(det).T
    bg = P.sum(axis=0)
    bg = (bg / bg.sum()).to_numpy()
    P = P.div(P.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    D = D.div(pd.Series(n_cells).replace(0, np.nan), axis=0).fillna(0)
    edges = np.arange(n_out + 1) * bin_size * rebin
    return dict(prof=P, bg=bg, detect=D, n_cells=n_cells, edges=edges)


def enrichment(res, pseudo_scale=0.05):
    pseudo = pseudo_scale / res["prof"].shape[1]
    return np.log2((res["prof"] + pseudo) / (res["bg"][None, :] + pseudo))


def call_regions(res, enr, min_log2=1.0, min_cell_frac=0.10, merge=True):
    rows, edges = [], res["edges"]
    for k in enr.index:
        hit = (enr.loc[k].to_numpy() >= min_log2) & (res["detect"].loc[k].to_numpy() >= min_cell_frac)
        idx = np.flatnonzero(hit)
        if len(idx) == 0:
            continue
        groups = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1) if merge else [[i] for i in idx]
        for g in groups:
            rows.append(dict(cluster=k, start=int(edges[g[0]]), end=int(edges[g[-1] + 1]),
                             n_bins=len(g),
                             max_log2_enrichment=float(enr.loc[k].iloc[g].max()),
                             mean_cell_frac=float(res["detect"].loc[k].iloc[g].mean()),
                             n_cells=res["n_cells"][k]))
    return pd.DataFrame(rows)


def plot_profile(res, enr, outfile, chrom="chr", cluster_names=None, vmax=2.0):
    edges = res["edges"]
    mid = (edges[:-1] + edges[1:]) / 2 / 1e6
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(14, 2.2 + 0.32 * len(enr)), sharex=True,
        gridspec_kw=dict(height_ratios=[1, max(2, 0.35 * len(enr))]))
    ax0.fill_between(mid, res["bg"], step="mid", color="0.4")
    ax0.set_ylabel("all cells\n(fraction)")
    ax0.set_title(f"{chrom}: Tn5 insertions ({(edges[1]-edges[0])/1e3:.0f} kb windows)")
    im = ax1.imshow(enr.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                    extent=[edges[0] / 1e6, edges[-1] / 1e6, len(enr), 0], interpolation="nearest")
    ax1.set_yticks(np.arange(len(enr)) + 0.5)
    ax1.set_yticklabels([str((cluster_names or {}).get(k, k)) for k in enr.index])
    ax1.set_xlabel("position (Mb)")
    ax1.set_ylabel("cluster")
    fig.colorbar(im, ax=[ax0, ax1], pad=0.01, label="log2(cluster / all cells)")
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(insertions, clusters, chrom="chr1", outdir=".", n_bins=None, bin_size=3000, rebin=10,
        min_log2=1.0, min_cell_frac=0.10, cluster_names=None, exclude_clusters=()):
    os.makedirs(outdir, exist_ok=True)
    clusters = clusters[~clusters.astype(str).isin([str(c) for c in exclude_clusters])]
    res = build_profile(insertions, clusters, n_bins, bin_size, rebin)
    enr = enrichment(res)
    plot_profile(res, enr, f"{outdir}/{chrom}_tn5_enrichment.png", chrom, cluster_names)
    table = call_regions(res, enr, min_log2, min_cell_frac)
    table.insert(0, "chrom", chrom)
    table.to_csv(f"{outdir}/{chrom}_enriched_regions.csv", index=False)
    return res, enr, table
