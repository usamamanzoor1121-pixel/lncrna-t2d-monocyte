#!/usr/bin/env python3
"""
pseudobulk.py
Utilities for pseudobulk aggregation from scRNA-seq data.
Used in bulk validation and co-expression analysis.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

def build_pseudobulk(adata, group_by=('sample_id', 'mono_subtype'),
                     layer='counts', min_cells=5):
    """
    Sum counts per (sample x group) pseudobulk.
    Returns: counts DataFrame (pseudobulks x genes), metadata DataFrame
    """
    counts_list, meta_list = [], []
    for keys, idx in adata.obs.groupby(list(group_by)).groups.items():
        if len(idx) < min_cells:
            continue
        X = adata[idx, :].layers.get(layer, adata[idx, :].X)
        if sp.issparse(X):
            s = np.array(X.sum(axis=0)).flatten()
        else:
            s = X.sum(axis=0)
        pb_id = '_'.join(str(k) for k in (keys if isinstance(keys, tuple) else [keys]))
        counts_list.append(pd.Series(s, index=adata.var_names, name=pb_id))
        meta_list.append(dict(zip(group_by, keys if isinstance(keys, tuple) else [keys])))
    counts_df = pd.DataFrame(counts_list).astype(int)
    meta_df   = pd.DataFrame(meta_list, index=counts_df.index)
    return counts_df, meta_df
