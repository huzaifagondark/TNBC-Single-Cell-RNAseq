import scanpy as sc
from pathlib import Path

BASE = Path(".")
RESULTS = BASE / "results"
FIGURES = BASE / "figures"
FIGURES.mkdir(exist_ok=True)

adata = sc.read_h5ad(RESULTS / "step1_raw_adata.h5ad")

# Mark mitochondrial genes
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

print(adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].describe())

# Save QC plots
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_qc_violin.png"
)

# Filtering thresholds
adata = adata[adata.obs["n_genes_by_counts"] > 200, :]
adata = adata[adata.obs["n_genes_by_counts"] < 7000, :]
adata = adata[adata.obs["pct_counts_mt"] < 10, :]

print("After QC:", adata)

adata.write(RESULTS / "step2_qc_adata.h5ad")
print("Saved:", RESULTS / "step2_qc_adata.h5ad")