import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

RESULTS = Path("results")
FIGURES = Path("figures")

adata = sc.read_h5ad(RESULTS / "step6_annotated_adata.h5ad")

# Cell type counts by patient
prop = (
    adata.obs
    .groupby(["patient_id", "cell_type"])
    .size()
    .reset_index(name="count")
)

total = prop.groupby("patient_id")["count"].transform("sum")
prop["percent"] = prop["count"] / total * 100

prop.to_csv(RESULTS / "cell_type_proportions_by_patient.csv", index=False)

# Plot
pivot = prop.pivot(index="patient_id", columns="cell_type", values="percent").fillna(0)

pivot.plot(kind="bar", stacked=True, figsize=(10, 6))
plt.ylabel("Percentage of cells")
plt.title("Cell type proportions by patient")
plt.tight_layout()
plt.savefig(FIGURES / "cell_type_proportions_by_patient.png", dpi=300)
plt.close()

print("Saved cell proportion table and plot.")