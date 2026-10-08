import pandas as pd
import scanpy as sc
from pathlib import Path

BASE = Path(".")
RAW = BASE / "data_raw"
CLEAN = BASE / "data_clean"
RESULTS = BASE / "results"

CLEAN.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

count_path = RAW / "GSE118389_counts_rsem.txt"

# Load count matrix: genes x cells
counts = pd.read_csv(count_path, sep="\t", index_col=0)

print("Counts shape genes x cells:", counts.shape)
print(counts.iloc[:5, :5])

# Build metadata from cell names
cell_ids = counts.columns.tolist()

def parse_cell_id(cell_id):
    parts = cell_id.split("_")
    patient_id = parts[0] if len(parts) > 0 else "unknown"
    plate_or_sample = parts[1] if len(parts) > 1 else "unknown"
    return {
        "cell_id": cell_id,
        "patient_id": patient_id,
        "sample_group": plate_or_sample
    }

metadata = pd.DataFrame([parse_cell_id(x) for x in cell_ids])
metadata.to_csv(CLEAN / "metadata_basic.csv", index=False)

print("Metadata shape:", metadata.shape)
print(metadata.head())

# Create AnnData
adata = sc.AnnData(counts.T)
adata.obs = metadata.set_index("cell_id")
adata.var_names_make_unique()

print("AnnData:", adata)

adata.write(RESULTS / "step1_raw_adata.h5ad")
print("Saved:", RESULTS / "step1_raw_adata.h5ad")