#!/usr/bin/env python3
"""
R_ANALYSIS Phase 0 — export data from the Python/anndata world into flat
files R can read directly, WITHOUT performing any statistical modeling here.
All actual inference (differential state, meta-analysis, pathway ranking,
deconvolution) happens in R on these exports, so the R results are a
genuinely independent computation, not a re-skin of the Python numbers.

Outputs (under scripts/r_analysis/data/):
  pseudobulk_counts_donor_subtype.csv   -- genes x (donor__subtype) raw counts
  pseudobulk_meta_donor_subtype.csv     -- donor, subtype, n_cells per pseudobulk sample
  gene_annotation.csv                   -- gene_symbol, ensembl_id, is_lncrna_gencode, mean_counts
  cell_level_key_genes.csv              -- per-cell: sample_id, mono_subtype, dpt_pseudotime,
                                            AC020656.1 and NEAT1 log_norm expr + raw counts
  bulk_counts.csv, bulk_meta.csv        -- GSE221521 counts/metadata (same parse as 04_bulk_validation.py)
  celltype_signature_matrix.csv         -- mean log_norm expr of marker genes per broad cell type,
                                            for R-side NNLS deconvolution
"""
import sys
import gzip
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

MONO_H5AD = "data/processed/GSE268210_monocytes_final.h5ad"
FULL_H5AD = "data/processed/GSE268210_phase1_full.h5ad"
GENCODE_LNC_CSV = "results/tables/lncrna_annotation_gencode_v32.csv"
OUTDIR = "scripts/r_analysis/data"

KEY_GENES = ["AC020656.1", "NEAT1", "MALAT1"]

# a small, standard PBMC marker panel (protein-coding only) for NNLS deconvolution
CELLTYPE_MARKERS = {
    "T_cell": ["CD3D", "CD3E", "CD3G", "IL7R", "CCR7"],
    "NK": ["NCAM1", "NKG7", "GNLY", "KLRF1", "KLRD1"],
    "Monocyte": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "VCAN"],
    "B_cell": ["CD79A", "CD79B", "MS4A1"],
    "DC": ["CD1C", "FCER1A", "LILRA4", "CLEC4C"],
    "Platelet": ["PPBP", "PF4", "GP1BA"],
}


def main():
    print("=== Loading monocyte h5ad ===", file=sys.stderr)
    mono = ad.read_h5ad(MONO_H5AD)
    print(f"  {mono.shape}", file=sys.stderr)

    # ---- gene annotation (from Phase 2 GENCODE output) ----
    gencode_lnc = pd.read_csv(GENCODE_LNC_CSV)
    lnc_set = set(gencode_lnc["gene_symbol"])
    gene_anno = mono.var.copy()
    gene_anno["gene_symbol"] = gene_anno.index
    gene_anno["is_lncrna_gencode"] = gene_anno["gene_symbol"].isin(lnc_set)
    gene_anno[["gene_symbol", "ensembl_id", "is_lncrna_gencode", "mean_counts", "n_cells"]].to_csv(
        f"{OUTDIR}/gene_annotation.csv", index=False)
    print(f"  Wrote gene_annotation.csv ({gene_anno.shape[0]} genes)", file=sys.stderr)

    # ---- donor x subtype pseudobulk RAW counts ----
    print("=== Building donor x subtype pseudobulk (raw counts) ===", file=sys.stderr)
    counts = mono.layers["counts"]
    donors = mono.obs["sample_id"].astype(str).values
    subtypes = mono.obs["mono_subtype"].astype(str).values
    genes = mono.var_names.tolist()

    pb_cols, pb_meta = [], []
    pb_matrix = []
    for donor in sorted(set(donors)):
        for st in ["CD14_Classical", "Intermediate", "CD16_NonClassical"]:
            mask = (donors == donor) & (subtypes == st)
            n = int(mask.sum())
            if n < 10:
                print(f"  SKIP {donor} x {st}: only {n} cells (<10)", file=sys.stderr)
                continue
            X = counts[mask, :]
            s = np.asarray(X.sum(axis=0)).flatten() if sp.issparse(X) else X.sum(axis=0)
            pb_matrix.append(s)
            col_name = f"{donor}__{st}"
            pb_cols.append(col_name)
            pb_meta.append({"pseudobulk_id": col_name, "donor": donor, "subtype": st, "n_cells": n,
                             "library_size": int(s.sum())})

    pb_df = pd.DataFrame(np.array(pb_matrix).T, index=genes, columns=pb_cols).astype(int)
    pb_df.to_csv(f"{OUTDIR}/pseudobulk_counts_donor_subtype.csv")
    pd.DataFrame(pb_meta).to_csv(f"{OUTDIR}/pseudobulk_meta_donor_subtype.csv", index=False)
    print(f"  Pseudobulk matrix: {pb_df.shape[0]} genes x {pb_df.shape[1]} donor x subtype samples", file=sys.stderr)
    print(f"  Samples per subtype:\n{pd.DataFrame(pb_meta)['subtype'].value_counts()}", file=sys.stderr)

    # ---- cell-level key gene table (for donor-level plots / sensitivity in R) ----
    print("=== Exporting cell-level key gene table ===", file=sys.stderr)
    rows = {"sample_id": donors, "mono_subtype": subtypes,
            "dpt_pseudotime": mono.obs["dpt_pseudotime"].values}
    for g in KEY_GENES:
        if g not in mono.var_names:
            continue
        expr = mono[:, g].layers["log_norm"]
        expr = np.asarray(expr.todense()).flatten() if sp.issparse(expr) else np.asarray(expr).flatten()
        raw = mono[:, g].layers["counts"]
        raw = np.asarray(raw.todense()).flatten() if sp.issparse(raw) else np.asarray(raw).flatten()
        rows[f"{g}_lognorm"] = expr
        rows[f"{g}_counts"] = raw
        rows[f"{g}_detected"] = (raw > 0).astype(int)
    cell_df = pd.DataFrame(rows)
    cell_df.to_csv(f"{OUTDIR}/cell_level_key_genes.csv", index=False)
    print(f"  Cell-level table: {cell_df.shape}", file=sys.stderr)

    del mono

    # ---- cell-type signature matrix from the FULL PBMC object (for NNLS deconvolution) ----
    print("=== Building cell-type signature matrix (full PBMC object) ===", file=sys.stderr)
    full = ad.read_h5ad(FULL_H5AD, backed="r")
    all_markers = sorted({g for gs in CELLTYPE_MARKERS.values() for g in gs if g in full.var_names})
    full_sub = full[:, all_markers].to_memory()
    sig_rows = {}
    for ct in CELLTYPE_MARKERS:
        mask = full_sub.obs["broad_celltype"] == ct
        if mask.sum() < 50:
            continue
        X = full_sub[mask, :].layers["log_norm"]
        m = np.asarray(X.mean(axis=0)).flatten() if sp.issparse(X) else X.mean(axis=0)
        sig_rows[ct] = m
    sig_df = pd.DataFrame(sig_rows, index=all_markers)
    sig_df.to_csv(f"{OUTDIR}/celltype_signature_matrix.csv")
    print(f"  Signature matrix: {sig_df.shape}", file=sys.stderr)

    print("\nDone. All exports in:", OUTDIR)


if __name__ == "__main__":
    main()
