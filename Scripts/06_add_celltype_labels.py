import scanpy as sc
from pathlib import Path

RESULTS = Path("results")

adata = sc.read_h5ad(RESULTS / "step5_markers_adata.h5ad")

# EDIT THIS after checking your dotplot and marker genes
cluster_to_celltype = {
    "0": "Unknown",
    "1": "Unknown",
    "2": "Unknown",
    "3": "Unknown",
    "4": "Unknown",
    "5": "Unknown",
    "6": "Unknown",
    "7": "Unknown",
    "8": "Unknown",
    "9": "Unknown"
}

adata.obs["cell_type"] = adata.obs["cluster"].map(cluster_to_celltype).fillna("Unknown")

sc.pl.umap(adata, color=["cell_type", "cluster"], save="_celltype_annotation.png")

adata.obs.to_csv(RESULTS / "cell_metadata_with_annotations.csv")
adata.write(RESULTS / "step6_annotated_adata.h5ad")

print(adata.obs[["cluster", "cell_type"]].value_counts())