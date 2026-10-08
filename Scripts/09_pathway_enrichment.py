import pandas as pd
import scanpy as sc
import gseapy as gp
from pathlib import Path
import re
RESULTS = Path("results")
INPUT_FILE = RESULTS / "step5_markers_adata.h5ad"
if not INPUT_FILE.exists():
    INPUT_FILE = RESULTS / "step4_clustered_adata.h5ad"
adata = sc.read_h5ad(INPUT_FILE)
GROUPBY = "cluster"
print("Using groupby:", GROUPBY)
print(adata.obs[GROUPBY].value_counts())
sc.tl.rank_genes_groups(
    adata,
    groupby=GROUPBY,
    method="wilcoxon",
    use_raw=True if adata.raw is not None else False
)
all_deg = []
groups = adata.uns["rank_genes_groups"]["names"].dtype.names
for group in groups:
    df = sc.get.rank_genes_groups_df(adata, group=group)
    df["group"] = group
    all_deg.append(df)
deg = pd.concat(all_deg, ignore_index=True)
deg.to_csv(RESULTS / "deg_by_cluster_full.csv", index=False)
print("DEG rows:", deg.shape)
# More relaxed filter for small dataset
deg_sig = deg[
    (deg["scores"] > 0) &
    (deg["logfoldchanges"] > 0)
].copy()
deg_sig.to_csv(RESULTS / "deg_by_cluster_positive.csv", index=False)
print("Positive DEG rows:", deg_sig.shape[0])
if deg_sig.empty:
    raise ValueError("No positive marker genes found. Check clustering/adata.")
def safe_name(name):
    return re.sub(r"[^\w\-]+", "_", str(name)).strip("_")
for group in deg_sig["group"].unique():
    genes = (
        deg_sig[deg_sig["group"] == group]
        .sort_values("scores", ascending=False)["names"]
        .dropna()
        .astype(str)
        .head(100)
        .tolist()
    )
    if len(genes) < 5:
        print(f"Skipping cluster {group}: fewer than 5 genes")
        continue
    print(f"Running enrichment for cluster {group}: {len(genes)} genes")
    try:
        gp.enrichr(
            gene_list=genes,
            gene_sets=[
                "GO_Biological_Process_2023",
                "KEGG_2021_Human",
                "MSigDB_Hallmark_2020"
            ],
            organism="Human",
            outdir=str(RESULTS / f"enrichment_cluster_{safe_name(group)}"),
            cutoff=0.1
        )
    except Exception as e:
        print(f"Failed cluster {group}: {e}")
print("Done.")