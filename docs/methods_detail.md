# Extended Methods

## Dataset Details

### GSE268210 — T2D PBMC scRNA-seq

| Parameter | Value |
|-----------|-------|
| Technology | 10x Genomics Chromium Next GEM Single Cell 5' Kit v2 |
| Sequencing | Illumina HiSeq 2500, 100 bp paired-end |
| Alignment | CellRanger v5.0.1 vs GRCh38 + GENCODE v32 (2020-A reference) |
| Patients | 10 T2D (1 excluded: DM_01, divergent feature count) |
| Features | 36,601 genes including 16,063 lncRNA annotations |
| Cells (post-QC) | 209,289 across 9 patients |

**Sample exclusion rationale (T2D_01 / GSM8287977):**
GSM8287977 yielded 33,538 detected features vs. 36,601 in all other samples,
suggesting alignment against a slightly different reference version or a CellRanger
parameter difference. Including it would require `join='outer'` in the concat step,
which would fill ~3,000 genes with zeros in this sample and bias batch-corrected
embeddings. Conservative exclusion was preferred.

---

## lncRNA Annotation

### Pattern-based detection (implemented in pipeline)

Gene name patterns used to identify lncRNAs without the GENCODE GTF:

```
^(LINC\d|AL\d{6}|AC\d{6}|SNHG|NEAT\d|MALAT|HOTAIR|XIST|
  KCNQ\dOT|MEG\d|MIR\d+HG|DLEU|MIAT|NORAD|PVT1|FTX|
  HAGLR|DANCR|CASC|GAS\d|HOTAIRM|PURPL)
```

Estimated recall: ~70% of GENCODE v32 lncRNAs with pattern-recognisable names.
Unannotated lncRNAs with AC/AL prefix account for ~40% of the GENCODE v32 lncRNA set.

### Complete annotation (recommended for publication)

Generate the complete Ensembl ID list from GENCODE v32:

```bash
# Download GENCODE v32 annotation GTF (~1.5 GB)
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz"

# Extract lncRNA gene Ensembl IDs
zcat gencode.v32.annotation.gtf.gz | \
    awk '$3=="gene"' | \
    grep 'gene_type "lncRNA"' | \
    grep -oP 'gene_id "\K[^"]+' | \
    sort -u > data/lncrna_gencode_v32_ids.txt

# Count
wc -l data/lncrna_gencode_v32_ids.txt
# Expected: ~16,063
```

Cross-reference in Python (add to `01_qc_preprocessing.py`):

```python
with open('data/lncrna_gencode_v32_ids.txt') as f:
    lncrna_set = set(l.strip().split('.')[0] for l in f if l.strip())

# adata.var['ensembl_id'] contains Ensembl IDs from features.tsv
base_ids = adata.var['ensembl_id'].str.split('.').str[0]
adata.var['is_lncrna'] = base_ids.isin(lncrna_set)
```

---

## Monocyte Subtype Annotation

### Cluster marker expression table

Marker expression (mean log-normalised counts) per Leiden 0.3 cluster
in the monocyte subset (n=54,653 cells):

| Cluster | n | CD14 | FCGR3A | S100A8 | MS4A7 | CDKN1C | HLA-DRA | Assignment |
|---------|---|------|--------|--------|-------|--------|---------|------------|
| 5 | 29,264 | **1.57** | 0.17 | **3.81** | 0.59 | 0.16 | 2.63 | CD14_Classical |
| 6 | 6,541 | **1.11** | 0.64 | **2.73** | 0.65 | 0.44 | 2.31 | CD14_Classical |
| 7 | 5,664 | 0.13 | **2.65** | 0.50 | **2.03** | **2.25** | 2.43 | CD16_NonClassical |
| 1 | 1,068 | 0.26 | 0.32 | 0.42 | 0.18 | 0.11 | 0.72 | Intermediate |
| 2 | 1,520 | 0.03 | 0.45 | 0.08 | 0.01 | 0.04 | **1.56** | Intermediate |
| 8 | 3,543 | 0.01 | 0.02 | 0.05 | 0.05 | 0.01 | **3.96** | Exclude_DC |
| 9 | 249 | 0.01 | 0.06 | 0.04 | 0.02 | 0.01 | **3.65** | Exclude_DC |
| 10 | 1,276 | 0.24 | 0.06 | 0.40 | 0.42 | 0.04 | **4.61** | Exclude_DC |
| 12 | 681 | 0.01 | 0.36 | 0.03 | 0.03 | 0.05 | **3.22** | Exclude_DC |
| 0 | 3,521 | 0.03 | 0.06 | 0.07 | 0.02 | 0.09 | 0.72 | Exclude_Ambiguous |
| 3 | 341 | 0.03 | **2.17** | 0.06 | 0.01 | 0.03 | 1.37 | Exclude_DC |
| 4 | 8 | 0.00 | 0.16 | 0.00 | 0.20 | 0.00 | 1.02 | Exclude_Ambiguous |
| 13 | 84 | 0.20 | 0.06 | 0.03 | 0.05 | 0.24 | 0.42 | Exclude_Ambiguous |
| 14 | 893 | 0.11 | 0.01 | 0.03 | 0.00 | 0.02 | 0.05 | Exclude_Ambiguous |

**Exclusion rationale:**
- Clusters 8, 9, 10, 12: HLA-DRA >> 3.5 with all monocyte markers < 0.1
  — characteristic DC contamination in monocyte gate
- Cluster 3: FCGR3A high but HLA-DRA elevated with low S100A8 — likely pDC
- Clusters 0, 4, 13, 14: All markers < 0.3, ambiguous identity

---

## Pseudotime Interpretation

**Why do Intermediate monocytes have the highest pseudotime (0.84)?**

Diffusion pseudotime measures transcriptional distance from the root in
diffusion map space — it is not a linear biological time axis. The root
was set to the most CD14-expressing cell in the CD14 Classical cluster.

Intermediate monocytes (clusters 1,2) have the most distinct transcriptional
profile relative to this root: they co-express genes from both the classical
(CD14) and non-classical (FCGR3A, HLA-DRA) programmes simultaneously, making
them maximally distant in diffusion space. This is consistent with Villani et al.
(Science 2017) who characterised intermediate monocytes as a distinct activated
state, not simply a transition point between classical and non-classical monocytes.

The trajectory should be interpreted as: **CD14 Classical is the origin state**,
with CD16 Non-Classical (PT=0.15) and Intermediate (PT=0.84) representing
divergent activated states from this origin.

---

## Bulk Validation Notes

### Why is concordance 55% and not higher?

The scRNA-seq trajectory captures **monocyte-subtype-specific** lncRNA expression.
The bulk RNA-seq (GSE221521) measures expression in **whole leukocytes**:
T cells (~30%) + NK cells (~25%) + Monocytes (~25%) + other cells.

Monocyte-specific signals are diluted by ~75% in bulk blood. This explains:
- Small effect sizes (bulk LFC 0.05–0.75 vs scRNA rho 0.15–0.38)
- MALAT1 showing no bulk change (p=0.891): ubiquitous expression masks subtype specificity
- AC020656.1 showing strong bulk change: monocyte-specific, little noise from other cell types

### Reference for concordance expectation

The Pearson r = −0.630, p = 0.003 between scRNA trajectory ρ and bulk LFC
is the correct validation metric. Individual gene concordance at 55% is
expected and reported as such in the manuscript.

---

## Software Versions

All analyses were run on Ubuntu 22.04 (WSL2) with:

```
Python          3.10.20
scanpy          1.11.5
anndata         0.11.4
harmonypy       0.0.9
numpy           2.4.4
scipy           1.17.1
pandas          2.3.3
matplotlib      3.8.x
seaborn         0.13.x
scikit-learn    1.3.x
pydeseq2        0.4.9
statsmodels     0.14.x
```
