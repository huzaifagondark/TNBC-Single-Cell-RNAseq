# Single-Cell Transcriptomics of Breast Cancer

**An exploratory, reproducible scRNA-seq portfolio analysis using Python and Scanpy**

**Author:** Huzaifa Salman  
**Project type:** Computational biology / cancer transcriptomics  
**Status:** Exploratory analysis completed; cell-type annotation and downstream comparisons require revision

## Overview

This project applies a single-cell RNA-sequencing (scRNA-seq) analysis workflow to a publicly available breast-cancer expression matrix containing **1,534 cells and 21,785 genes**. The aim is to explore cellular heterogeneity and tumor-microenvironment transcriptional programs while documenting the distinction between computational outputs and biologically validated conclusions.

The workflow covers quality-control visualization, library-size normalization and log transformation, highly variable gene selection, principal component analysis (PCA), neighborhood graph construction, UMAP, Leiden clustering, canonical marker analysis, cluster-level marker ranking, and exploratory gene-signature scoring.

> **Important:** This is an exploratory portfolio project, not a validated clinical analysis or an exact reproduction of the reference study. The current `cell_type` metadata column contains only `Unknown`; therefore, existing cell-type proportion and `Unknown vs. rest` differential-expression plots are **not biologically interpretable**.

## Key findings

- **16 Leiden clusters** (clusters 0–15) were identified.
- Canonical marker patterns suggest **epithelial-like, T-cell, B-cell, NK/cytotoxic, myeloid/macrophage, fibroblast/stromal, and endothelial-like** populations. These are **provisional marker-based interpretations**, not finalized cell-type annotations.
- Patient- and sample-group-colored UMAPs show both shared transcriptional structure and some sample-associated regions, indicating that patient/batch effects deserve further evaluation.
- Exploratory gene-signature maps show heterogeneous expression patterns related to **T-cell exhaustion, regulatory T cells, macrophage activity, TNF/NF-κB inflammation, fibroblast extracellular matrix, and hypoxia**. These scores do not establish pathway activation or clinical associations.

## Workflow

```text
Processed expression matrix (21,785 genes × 1,534 cells)
        ↓
AnnData construction and metadata organization
        ↓
Quality-control metrics and visualization
        ↓
Library-size normalization + log1p
        ↓
Highly variable gene selection
        ↓
PCA → k-nearest-neighbor graph
        ↓
UMAP → Leiden clustering (16 clusters)
        ↓
Canonical marker visualization and cluster-level marker ranking
        ↓
Exploratory tumor-microenvironment signature scoring
        ↓
Audit of annotation and downstream-analysis limitations
```

## Methods and tools

| Analysis stage | Approach |
| --- | --- |
| Data handling | Python, Pandas, AnnData |
| Quality control | Scanpy QC metrics and distribution plots |
| Normalization | Library-size normalization, `log1p` |
| Feature selection | Highly variable genes |
| Dimensionality reduction | PCA |
| Neighborhood structure | k-nearest-neighbor graph |
| Visualization | UMAP |
| Clustering | Leiden algorithm |
| Marker identification | Canonical markers and Wilcoxon-based cluster marker ranking |
| Functional exploration | Relative gene-signature scoring |

## Biological marker panels

| Putative population | Representative markers |
| --- | --- |
| Epithelial-like | `EPCAM`, `KRT8`, `KRT18`, `KRT19` |
| T cells | `CD3D`, `CD3E` |
| B cells | `MS4A1`, `CD79A`, `CD79B` |
| NK/cytotoxic lymphocytes | `NKG7`, `GNLY`, `KLRD1` |
| Myeloid/macrophage-like | `LYZ`, `LST1`, `CD68`, `CST3` |
| Fibroblast/stromal | `COL1A1`, `COL1A2`, `DCN`, `LUM` |
| Endothelial-like | `PECAM1`, `VWF`, `KDR` |

Marker expression supports these hypotheses, but mixed/ambiguous clusters require additional validation before labels are finalized.

## Figures and outputs

The accompanying portfolio report documents:

1. QC distributions for detected genes, total counts, and mitochondrial percentage
2. PCA explained-variance ranking
3. UMAP colored by Leiden cluster
4. UMAP colored by patient ID and sample group
5. Canonical marker-gene feature plots
6. Cluster-level marker dot plot
7. Ranked marker genes for clusters 0–15
8. Exploratory immune, stromal, inflammatory, and hypoxia signature maps

**Note:** The report also includes cell-type annotation, cell-type proportions, and `Unknown vs. rest` differential-expression figures specifically to document the current annotation issue. These should **not** be interpreted as biological results.

## Repository organization

```text
.
├── Figures/                  # Visualizations
├── Results/                  # Analysis outputs
├── Scripts/                  # Analysis scripts
├── data_clean/               # Locally prepared data, if distributable
└── README.md
```

A detailed PDF portfolio report may also be included in the repository. Directory contents and filenames may vary. For reproducibility, review the scripts in `Scripts/` and document the original dataset accession, file provenance, environment/package versions, and execution order before treating this as a fully reproducible pipeline.

## Limitations

- **Incomplete annotation:** All entries in `cell_type` are currently `Unknown`, so cell-type composition and cell-type-level differential expression need to be rerun after validated annotation.
- **Mitochondrial QC:** The computed mitochondrial percentage is uniformly zero under the current `MT-` gene-name rule; it is **not** evidence of zero mitochondrial RNA or absence of stressed cells.
- **Patient/sample effects:** Some clusters show patient/sample enrichment, warranting investigation of technical and biological variation.
- **Exploratory signatures:** Gene-set scores are relative and do not establish pathway activation.
- **Reference study:** This project follows core single-cell analytical ideas but does not reproduce the original study's full cohort or advanced analyses.

## Next steps

- [ ] Finalize cluster-to-cell-type labels using marker evidence; retain `Unresolved` for ambiguous clusters
- [ ] Save corrected AnnData metadata and regenerate cell-type UMAPs
- [ ] Recalculate cell-type proportions by patient and sample
- [ ] Perform valid cell-type/condition comparisons with an appropriate statistical design
- [ ] Explore enrichment with multiple-testing correction
- [ ] Evaluate patient/sample effects and integration methods where justified
- [ ] Document data accession, software versions, and step-by-step reproduction instructions

## Project report

See the **Single-Cell Transcriptomics Portfolio Report** PDF in this project for figures, interpretation, and a detailed discussion of current limitations.

## Attribution

This is an independent educational/portfolio analysis of publicly available breast-cancer expression data. The source dataset and reference study should be cited explicitly once their identifiers and bibliographic details are confirmed from the project files.

---

**Huzaifa Salman**  
*Bioinformatics | Single-Cell Transcriptomics | Cancer Biology*
