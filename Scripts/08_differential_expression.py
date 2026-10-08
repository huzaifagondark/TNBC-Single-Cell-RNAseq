import scanpy as sc
from pathlib import Path

RESULTS = Path("results")

adata = sc.read_h5ad(RESULTS / "step6_annotated_adata.h5ad")

# DEG by cell type
sc.tl.rank_genes_groups(
    adata,
    groupby="cell_type",
    method="wilcoxon",
    use_raw=True
)

sc.pl.rank_genes_groups(adata, n_genes=15, sharey=False, save="_deg_celltypes.png")

deg_df = sc.get.rank_genes_groups_df(adata, group=None)
deg_df.to_csv(RESULTS / "deg_by_cell_type.csv", index=False)

print("Saved:", RESULTS / "deg_by_cell_type.csv")