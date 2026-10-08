import scanpy as sc
from pathlib import Path

RESULTS = Path("results")

adata = sc.read_h5ad(RESULTS / "step4_clustered_adata.h5ad")

marker_genes = {
    "Malignant/Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T cells": ["CD3D", "CD3E", "TRAC"],
    "B cells": ["MS4A1", "CD79A", "CD79B"],
    "NK cells": ["NKG7", "GNLY", "KLRD1"],
    "Myeloid/Macrophage": ["LYZ", "LST1", "CD68", "CST3"],
    "Fibroblast/Stromal": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "KDR"]
}

# Use raw gene list if available
gene_source = adata.raw.var_names if adata.raw is not None else adata.var_names

# Filter missing genes from each marker group
filtered_marker_genes = {}
missing_genes = []

for cell_type, genes in marker_genes.items():
    present = [g for g in genes if g in gene_source]
    missing = [g for g in genes if g not in gene_source]
    if present:
        filtered_marker_genes[cell_type] = present
    missing_genes.extend(missing)

print("Missing marker genes removed:", missing_genes)
print("Markers used:", filtered_marker_genes)

# Make one flat marker list for UMAP
available_markers = []
for genes in filtered_marker_genes.values():
    available_markers.extend(genes)

# Plot marker genes on UMAP
sc.pl.umap(
    adata,
    color=available_markers,
    use_raw=True if adata.raw is not None else False,
    save="_marker_genes.png"
)

# Dotplot only with genes that exist
sc.pl.dotplot(
    adata,
    filtered_marker_genes,
    groupby="cluster",
    use_raw=True if adata.raw is not None else False,
    standard_scale="var",
    save="_marker_dotplot.png"
)

# Cluster marker genes
sc.tl.rank_genes_groups(
    adata,
    groupby="cluster",
    method="wilcoxon",
    use_raw=True if adata.raw is not None else False
)

sc.pl.rank_genes_groups(
    adata,
    n_genes=10,
    sharey=False,
    save="_cluster_markers.png"
)

markers_df = sc.get.rank_genes_groups_df(adata, group=None)
markers_df.to_csv(RESULTS / "cluster_marker_genes.csv", index=False)

adata.write(RESULTS / "step5_markers_adata.h5ad")

print("Done. Marker annotation files saved.")