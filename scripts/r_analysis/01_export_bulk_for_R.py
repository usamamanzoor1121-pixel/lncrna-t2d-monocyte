#!/usr/bin/env python3
"""
R_ANALYSIS Phase 0b — export GSE221521 bulk counts + metadata to flat CSV for
an independent R-based reanalysis (limma/edgeR), reusing the same parsing
logic already validated in scripts/analysis/04_bulk_validation.py (this is
data extraction only, not statistical modeling -- the actual DE testing is
left to R).
"""
import gzip
import numpy as np
import pandas as pd
from pathlib import Path

DATA_DIR = Path("data/bulk_validation/GSE221521_data")
OUTDIR = Path("scripts/r_analysis/data")


def load_counts():
    count_file = DATA_DIR / "GSE221521_gene_expression.xls.gz"
    full_df = pd.read_csv(count_file, index_col=0, sep="\t", compression="gzip")
    if "gene_name" in full_df.columns:
        full_df.index = full_df["gene_name"].values
        full_df = full_df.drop(columns=["gene_name"])
    anno_cols = [c for c in full_df.columns
                 if c in ["description", "gene_type", "locus", "chr", "start", "end", "strand", "length", "gene_id"]]
    full_df = full_df.drop(columns=anno_cols)
    col_names = [str(c).lower() for c in full_df.columns]
    count_cols = [c for c, cn in zip(full_df.columns, col_names) if "count" in cn]
    counts_df = full_df[count_cols].copy() if count_cols else full_df.copy()
    counts_df.columns = [str(c).replace("_count", "").replace("_FPKM", "") for c in counts_df.columns]
    counts_df = counts_df[~counts_df.index.duplicated(keep="first")]
    return counts_df


def parse_series_matrix_detailed(series_file, count_sample_ids):
    with gzip.open(series_file, "rt") as f:
        lines = f.readlines()
    titles = []
    for line in lines:
        line = line.strip()
        if line.startswith("!Sample_title"):
            titles = [t.strip('"') for t in line.split("\t")[1:]]
            break
    sample_to_condition = {}
    for title in titles:
        parts = title.strip().split()
        sid = parts[-1] if parts else ""
        t_low = title.lower()
        if " dm " in t_low or t_low.startswith("dm") or ", dm " in t_low:
            cond = "T2D"
        elif " dr " in t_low or t_low.startswith("dr") or ", dr " in t_low:
            cond = "PreDM"
        elif "control" in t_low:
            cond = "Control"
        else:
            cond = "Unknown"
        sample_to_condition[sid] = cond
    conditions = [sample_to_condition.get(c, "Unknown") for c in count_sample_ids]
    return pd.DataFrame({"condition": conditions}, index=count_sample_ids)


def main():
    counts_df = load_counts()
    series_file = DATA_DIR / "GSE221521_series_matrix.txt.gz"
    meta_df = parse_series_matrix_detailed(series_file, list(counts_df.columns))

    keep = meta_df["condition"] != "Unknown"
    counts_df = counts_df.loc[:, keep.values]
    meta_df = meta_df.loc[keep]

    counts_df.to_csv(OUTDIR / "bulk_counts.csv")
    meta_df.to_csv(OUTDIR / "bulk_meta.csv")
    print(f"Bulk counts: {counts_df.shape}  |  Conditions: {meta_df['condition'].value_counts().to_dict()}")
    print(f"Written to {OUTDIR}/bulk_counts.csv, bulk_meta.csv")


if __name__ == "__main__":
    main()
