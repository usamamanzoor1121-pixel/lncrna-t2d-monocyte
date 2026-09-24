#!/usr/bin/env python3
"""
05_paga_topology.py
Phase 5: Branch-aware trajectory topology (PAGA).

Diffusion pseudotime imposes a single linear ordering on the data by
construction, which can be misleading when the underlying structure is not
linear. PAGA (Partition-based Graph Abstraction) computes a connectivity
strength between each pair of monocyte subtypes on the same Harmony-corrected
embedding, making any branching structure explicit rather than assumed.

Usage:
  python scripts/analysis/05_paga_topology.py \
      --mono_h5ad data/processed/GSE268210_monocytes_final.h5ad \
      --outdir results/tables --figdir results/figures
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 1


def main():
    parser = argparse.ArgumentParser(description="PAGA connectivity between monocyte subtypes")
    parser.add_argument("--mono_h5ad", required=True)
    parser.add_argument("--outdir", default="results/tables")
    parser.add_argument("--figdir", default="results/figures")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.mono_h5ad)
    use_rep = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"
    print(f"Using representation: {use_rep}  shape={adata.obsm[use_rep].shape}")

    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15, random_state=42)
    adata.obs["mono_subtype"] = adata.obs["mono_subtype"].astype("category")
    sc.tl.paga(adata, groups="mono_subtype")

    conn = pd.DataFrame(
        adata.uns["paga"]["connectivities"].toarray(),
        index=adata.obs["mono_subtype"].cat.categories,
        columns=adata.obs["mono_subtype"].cat.categories,
    )
    conn.to_csv(f"{args.outdir}/paga_connectivity.csv")
    print("\nPAGA connectivity matrix (subtype x subtype):")
    print(conn.to_string())

    cats = list(adata.obs["mono_subtype"].cat.categories)
    edges = []
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            edges.append((cats[i], cats[j], conn.iloc[i, j]))
    edges.sort(key=lambda x: -x[2])
    print("\nPairwise connectivity strengths (higher = more direct transcriptional continuity):")
    for a, b, w in edges:
        print(f"  {a:<20} <-> {b:<20}  connectivity = {w:.4f}")

    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.paga(adata, color="mono_subtype", ax=ax, show=False, node_size_scale=2,
               fontsize=10, threshold=0.0, edge_width_scale=1.5)
    fig.suptitle("PAGA Connectivity -- Monocyte Subtypes (T2D)", fontsize=12, fontweight="bold")
    fig.savefig(f"{args.figdir}/Fig_paga_topology.png", dpi=150, bbox_inches="tight", facecolor="white")
    print(f"\nFigure saved: {args.figdir}/Fig_paga_topology.png")


if __name__ == "__main__":
    main()
