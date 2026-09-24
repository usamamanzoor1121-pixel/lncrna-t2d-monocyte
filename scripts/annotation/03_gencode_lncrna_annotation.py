#!/usr/bin/env python3
"""
03_gencode_lncrna_annotation.py
Phase 3: GENCODE v32 lncRNA annotation.

lncRNAs are identified by cross-referencing each gene's Ensembl ID (from the
10x feature reference) against the complete GENCODE v32 annotation GTF,
retaining genes with gene_type == "lncRNA". This is the authoritative
annotation used throughout the rest of the pipeline.

For reference, scripts/utils/lncrna_annotation.py also provides a fast
gene-symbol pattern-matching heuristic (regex on LINC*/AC*/AL*/SNHG*/NEAT*
prefixes, etc.). Benchmarked against the GENCODE cross-reference below, that
heuristic achieves 71.5% recall and 98.2% precision -- useful as a quick
first-pass filter, but the GENCODE cross-reference is what all downstream
lncRNA-specific analyses in this repository actually use.

Inputs:
  data/reference/gencode.v32.annotation.gtf.gz   (GENCODE v32 GTF; download command below)
  data/processed/GSE268210_monocytes_final.h5ad  (monocyte object, has var['ensembl_id'])

Outputs (results/tables/):
  gencode_v32_lncrna_ids.txt         -- all GENCODE v32 lncRNA Ensembl gene IDs (versionless)
  lncrna_annotation_gencode_v32.csv  -- per-gene table: gene_id, gene_name, gene_type, expression summary
  lncrna_annotation_pattern_vs_gencode.csv -- full pattern-heuristic vs. GENCODE comparison

Download the reference (once):
  wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz" \
       -O data/reference/gencode.v32.annotation.gtf.gz

Usage:
  python scripts/annotation/03_gencode_lncrna_annotation.py \
      --gtf data/reference/gencode.v32.annotation.gtf.gz \
      --h5ad data/processed/GSE268210_monocytes_final.h5ad \
      --outdir results/tables
"""
import argparse
import gzip
import re
import sys
import numpy as np
import pandas as pd
import anndata as ad


def parse_gencode_genes(gtf_path):
    """Parse gene-level GTF lines -> DataFrame[gene_id, gene_name, gene_type]."""
    rows = []
    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene_name "([^"]+)"')
    gene_type_re = re.compile(r'gene_type "([^"]+)"')
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            gid = gene_id_re.search(attrs)
            gname = gene_name_re.search(attrs)
            gtype = gene_type_re.search(attrs)
            rows.append({
                "gene_id_versioned": gid.group(1) if gid else None,
                "gene_name": gname.group(1) if gname else None,
                "gene_type": gtype.group(1) if gtype else None,
                "chrom": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
            })
    df = pd.DataFrame(rows)
    df["gene_id"] = df["gene_id_versioned"].str.split(".").str[0]
    return df


def main():
    parser = argparse.ArgumentParser(description="GENCODE v32 lncRNA annotation")
    parser.add_argument("--gtf", required=True, help="Path to gencode.v32.annotation.gtf.gz")
    parser.add_argument("--h5ad", required=True, help="Path to the monocyte h5ad (has var['ensembl_id'])")
    parser.add_argument("--outdir", default="results/tables")
    args = parser.parse_args()

    print("Parsing GENCODE v32 GTF...", file=sys.stderr)
    gtf_genes = parse_gencode_genes(args.gtf)
    print(f"  Total genes in GTF: {len(gtf_genes):,}", file=sys.stderr)

    lncrna_genes = gtf_genes[gtf_genes["gene_type"] == "lncRNA"].copy()
    print(f"  GENCODE v32 lncRNA genes: {len(lncrna_genes):,}", file=sys.stderr)

    lncrna_id_set = set(lncrna_genes["gene_id"])
    with open(f"{args.outdir}/gencode_v32_lncrna_ids.txt", "w") as f:
        for gid in sorted(lncrna_id_set):
            f.write(gid + "\n")

    print("Loading monocyte h5ad (var only, backed mode)...", file=sys.stderr)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    var = adata.var.copy()
    var["gene_symbol"] = var.index
    var["ensembl_base"] = var["ensembl_id"].astype(str).str.split(".").str[0]

    gtype_map = dict(zip(gtf_genes["gene_id"], gtf_genes["gene_type"]))
    gname_map = dict(zip(gtf_genes["gene_id"], gtf_genes["gene_name"]))
    var["gencode_gene_type"] = var["ensembl_base"].map(gtype_map)
    var["gencode_gene_name"] = var["ensembl_base"].map(gname_map)
    var["is_lncrna_gencode"] = var["gencode_gene_type"] == "lncRNA"
    var["found_in_gencode_v32"] = var["ensembl_base"].isin(gtf_genes["gene_id"])

    pattern_n = int(var["is_lncrna"].sum()) if "is_lncrna" in var.columns else 0
    gencode_n = int(var["is_lncrna_gencode"].sum())
    both = int((var.get("is_lncrna", False) & var["is_lncrna_gencode"]).sum())

    print("\n=== GENCODE v32 lncRNA annotation summary ===")
    print(f"  Dataset genes total                    : {var.shape[0]:,}")
    print(f"  GENCODE-v32-confirmed lncRNA count     : {gencode_n:,}")
    if pattern_n > 0:
        print(f"  Pattern-heuristic recall vs. GENCODE   : {both/gencode_n*100:.1f}%")
        print(f"  Pattern-heuristic precision             : {both/pattern_n*100:.1f}%")

    out_cols = ["gene_symbol", "ensembl_id", "ensembl_base", "is_lncrna",
                "found_in_gencode_v32", "gencode_gene_type", "gencode_gene_name",
                "is_lncrna_gencode", "mean_counts", "n_cells", "pct_dropout_by_counts"]
    out_cols = [c for c in out_cols if c in var.columns]
    var[out_cols].to_csv(f"{args.outdir}/lncrna_annotation_pattern_vs_gencode.csv", index=False)

    gencode_lnc = var[var["is_lncrna_gencode"]].copy()
    keep = ["gene_symbol", "ensembl_id", "gencode_gene_type", "mean_counts", "n_cells", "pct_dropout_by_counts"]
    keep = [c for c in keep if c in gencode_lnc.columns]
    gencode_lnc[keep].rename(columns={"gencode_gene_type": "gene_type"}).assign(
        annotation_source="GENCODE_v32_ensembl_crossref"
    ).to_csv(f"{args.outdir}/lncrna_annotation_gencode_v32.csv", index=False)

    print(f"\nOutputs written to {args.outdir}/")


if __name__ == "__main__":
    main()
