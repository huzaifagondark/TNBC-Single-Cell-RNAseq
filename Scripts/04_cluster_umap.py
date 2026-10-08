import scanpy as sc
from pathlib import Path

RESULTS = Path("results")

adata = sc.read_h5ad(RESULTS / "step3_normalized_hvg_adata.h5ad")

# Scale data
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver="arpack")

sc.pl.pca_variance_ratio(adata, n_pcs=40, log=True, save="_pca_variance.png")

# Neighbors and UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)

# Leiden clustering
sc.tl.leiden(adata, resolution=0.5, key_added="cluster")

# Plots
sc.pl.umap(adata, color=["cluster"], save="_clusters.png")
sc.pl.umap(adata, color=["patient_id"], save="_patients.png")
sc.pl.umap(adata, color=["sample_group"], save="_sample_group.png")

adata.write(RESULTS / "step4_clustered_adata.h5ad")
print("Saved:", RESULTS / "step4_clustered_adata.h5ad")