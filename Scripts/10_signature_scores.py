import scanpy as sc
from pathlib import Path

RESULTS = Path("results")

adata = sc.read_h5ad(RESULTS / "step6_annotated_adata.h5ad")

signatures = {
    "T_cell_exhaustion": ["PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2"],
    "Treg_signature": ["FOXP3", "IL2RA", "CTLA4"],
    "Macrophage_pro_tumor": ["CCL2", "SPP1", "MMP9", "MGP"],
    "Inflammation_TNF_NFKB": ["NFKB1", "RELA", "TNF", "CXCL8", "CCL20"],
    "Fibroblast_ECM": ["COL1A1", "COL1A2", "FN1", "DCN", "LUM"],
    "Hypoxia": ["VEGFA", "CA9", "LDHA", "ENO1", "SLC2A1"]
}

for name, genes in signatures.items():
    genes_present = [g for g in genes if g in adata.raw.var_names]
    if len(genes_present) >= 2:
        sc.tl.score_genes(adata, gene_list=genes_present, score_name=name, use_raw=True)
        print(name, genes_present)
    else:
        print("Skipped", name, "not enough genes present")

sc.pl.umap(
    adata,
    color=[s for s in signatures.keys() if s in adata.obs.columns],
    save="_signature_scores.png"
)

adata.obs.to_csv(RESULTS / "metadata_with_signature_scores.csv")
adata.write(RESULTS / "step10_signature_scores_adata.h5ad")