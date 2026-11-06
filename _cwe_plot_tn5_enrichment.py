#!/usr/bin/env python3
"""
Plot Tn5 insertion profiles across each chromosome, per ATAC cluster, to see
whether any regions are enriched in particular cell types.

Assumed inputs (edit `positions_by_chrom` below if your dictionary differs):
    insertions : dict  barcode -> dict  chrom -> 1-D array of insertion positions (bp)
    clusters   : pandas Series  barcode -> cluster label   (e.g. a.obs["atac_leiden"])

Method (mirrors what the per-chromosome entropy summarizes):
  * Each chromosome is split into fixed bins (default 100 kb).
  * Every cell contributes a profile that sums to 1 over the chromosome, so deeply
    sequenced cells do not dominate. Cell profiles are summed per cluster and
    renormalised -> cluster profile p_c(bin).
  * Background = profile of all cells pooled (equal weight per cell).
  * Enrichment = log2((p_c + pseudo) / (p_bg + pseudo)).
  * A bin is called enriched only if log2 enrichment >= --min-log2 AND at least
    --min-cell-frac of the cluster's cells have an insertion in that bin.
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------- input adapter
def positions_by_chrom(entry):
    """Return {chrom: np.ndarray of positions} for one cell's entry.
    EDIT THIS if your per-barcode value is not a dict of chrom -> positions."""
    if isinstance(entry, dict):
        return {c: np.asarray(v) for c, v in entry.items()}
    if hasattr(entry, "columns"):                      # DataFrame with chrom / pos columns
        pos_col = "pos" if "pos" in entry.columns else "start"
        return {c: g[pos_col].to_numpy() for c, g in entry.groupby("chrom")}
    raise TypeError("Unsupported insertion format; edit positions_by_chrom().")


# ---------------------------------------------------------------- core
def build_profiles(insertions, clusters, chrom_sizes, bin_size=100_000, chroms=None):
    """Return {chrom: dict(prof, bg, detect, n_cells, edges)} with
    prof/detect as DataFrames (clusters x bins)."""
    clusters = clusters[clusters.index.isin(insertions.keys())]
    labels = sorted(clusters.unique(), key=lambda x: str(x))
    chroms = chroms or list(chrom_sizes)
    out = {}
    for chrom in chroms:
        n_bins = int(np.ceil(chrom_sizes[chrom] / bin_size))
        edges = np.arange(n_bins + 1) * bin_size
        prof = {k: np.zeros(n_bins) for k in labels}       # sum of per-cell fractions
        det = {k: np.zeros(n_bins) for k in labels}        # cells with >=1 insertion
        n_cells = {k: 0 for k in labels}
        for bc, k in clusters.items():
            pos = positions_by_chrom(insertions[bc]).get(chrom)
            if pos is None or len(pos) == 0:
                continue
            counts, _ = np.histogram(pos, bins=edges)
            prof[k] += counts / counts.sum()               # each cell weighs 1
            det[k] += counts > 0
            n_cells[k] += 1
        P = pd.DataFrame(prof).T
        D = pd.DataFrame(det).T
        bg = P.sum(axis=0)
        bg = bg / bg.sum() if bg.sum() > 0 else bg
        P = P.div(P.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        D = D.div(pd.Series(n_cells).replace(0, np.nan), axis=0).fillna(0)
        out[chrom] = dict(prof=P, bg=bg.to_numpy(), detect=D, n_cells=n_cells, edges=edges)
    return out


def enrichment(res, pseudo_scale=0.05):
    pseudo = pseudo_scale / res["prof"].shape[1]
    return np.log2((res["prof"] + pseudo) / (res["bg"][None, :] + pseudo))


def call_regions(res, enr, min_log2=1.0, min_cell_frac=0.10):
    rows = []
    edges = res["edges"]
    for k in enr.index:
        hit = (enr.loc[k].to_numpy() >= min_log2) & (res["detect"].loc[k].to_numpy() >= min_cell_frac)
        for i in np.flatnonzero(hit):
            rows.append(dict(cluster=k, start=int(edges[i]), end=int(edges[i + 1]),
                             log2_enrichment=float(enr.loc[k].iloc[i]),
                             cell_frac=float(res["detect"].loc[k].iloc[i]),
                             n_cells=res["n_cells"][k]))
    return rows


def plot_chrom(chrom, res, enr, outfile, cluster_names=None, vmax=2.0):
    edges = res["edges"]
    mid_mb = (edges[:-1] + edges[1:]) / 2 / 1e6
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=(14, 2.2 + 0.32 * len(enr)), sharex=True,
        gridspec_kw=dict(height_ratios=[1, max(2, 0.35 * len(enr))]))
    ax0.fill_between(mid_mb, res["bg"], step="mid", color="0.4")
    ax0.set_ylabel("all cells\n(fraction)")
    ax0.set_title(f"{chrom}: Tn5 insertion profile ({(edges[1]-edges[0])/1e3:.0f} kb bins)")
    im = ax1.imshow(enr.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                    extent=[edges[0] / 1e6, edges[-1] / 1e6, len(enr), 0], interpolation="nearest")
    names = [str((cluster_names or {}).get(k, k)) for k in enr.index]
    ax1.set_yticks(np.arange(len(enr)) + 0.5)
    ax1.set_yticklabels(names)
    ax1.set_xlabel("position (Mb)")
    ax1.set_ylabel("cluster")
    fig.colorbar(im, ax=[ax0, ax1], pad=0.01, label="log2(cluster / all cells)")
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(insertions, clusters, chrom_sizes, outdir=".", bin_size=100_000, chroms=None,
        min_log2=1.0, min_cell_frac=0.10, cluster_names=None, exclude_clusters=()):
    import os
    os.makedirs(outdir, exist_ok=True)
    clusters = clusters[~clusters.astype(str).isin([str(c) for c in exclude_clusters])]
    profiles = build_profiles(insertions, clusters, chrom_sizes, bin_size, chroms)
    all_rows = []
    for chrom, res in profiles.items():
        enr = enrichment(res)
        plot_chrom(chrom, res, enr, f"{outdir}/{chrom}_tn5_enrichment.png", cluster_names)
        for r in call_regions(res, enr, min_log2, min_cell_frac):
            all_rows.append(dict(chrom=chrom, **r))
    table = pd.DataFrame(all_rows)
    table.to_csv(f"{outdir}/enriched_regions.csv", index=False)
    return profiles, table


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run from Python instead: see run().")
    ap.parse_args()
    print("Import and call run(insertions, clusters, chrom_sizes, outdir=...) from your session.")
