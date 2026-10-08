import scanpy as sc
from pathlib import Path

RESULTS = Path("results")

adata = sc.read_h5ad(RESULTS / "step2_qc_adata.h5ad")

# Save raw counts before normalization
adata.layers["counts"] = adata.X.copy()

# Normalize each cell to same total expression
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform
sc.pp.log1p(adata)

# Save normalized data
adata.raw = adata

# Highly variable genes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat"
)

print("Highly variable genes:", adata.var["highly_variable"].sum())

# Keep only HVGs for PCA/clustering
adata = adata[:, adata.var["highly_variable"]].copy()

adata.write(RESULTS / "step3_normalized_hvg_adata.h5ad")
print("Saved:", RESULTS / "step3_normalized_hvg_adata.h5ad")