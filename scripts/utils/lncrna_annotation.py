#!/usr/bin/env python3
"""
lncrna_annotation.py
Utilities for lncRNA gene annotation.
- Pattern-based detection from gene names (~70% recall)
- GENCODE v32 Ensembl ID cross-reference (complete, requires GTF)
"""
import re
import pandas as pd

LNCRNA_RE = re.compile(
    r'^(LINC\d|AL\d{6}|AC\d{6}|SNHG|NEAT\d|MALAT|HOTAIR|XIST|'
    r'KCNQ\dOT|MEG\d|MIR\d+HG|DLEU|MIAT|NORAD|PVT1|FTX|'
    r'HAGLR|DANCR|CASC|GAS\d|HOTAIRM|PURPL)',
    re.IGNORECASE
)

def annotate_by_name(gene_names):
    """Pattern-based lncRNA detection. Returns boolean Series."""
    return pd.Series(gene_names).map(lambda x: bool(LNCRNA_RE.match(str(x)))).values

def annotate_by_gencode(gene_names, ensembl_ids, lncrna_ids_file):
    """
    Cross-reference against GENCODE v32 lncRNA Ensembl IDs.
    More accurate than name patterns.
    Generate lncrna_ids_file: see docs/methods_detail.md
    """
    with open(lncrna_ids_file) as f:
        lncrna_set = set(l.strip().split('.')[0] for l in f if l.strip())
    base_ids = pd.Series(ensembl_ids).str.split('.').str[0]
    return base_ids.isin(lncrna_set).values
