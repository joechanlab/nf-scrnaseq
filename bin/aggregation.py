import argparse
import scanpy as sc
import numpy as np
from utils import get_basename_without_extension

parser = argparse.ArgumentParser(description="wrapper for concatenating the samples.")
parser.add_argument(
    "inputs", nargs="*", help="Paths to directory containing files to be concatenated."
)
parser.add_argument("--output", required=True, help="The output h5ad file path.")

args = parser.parse_args()

adata_list = []

for idx, h5_file_path in enumerate(args.inputs, start=1):
    print(h5_file_path)
    adata = sc.read(h5_file_path)
    adata.obs["sample_name"] = get_basename_without_extension(h5_file_path).rsplit(
        "_", 1
    )[0]
    adata.obs.index = adata.obs["sample_name"] + "_" + adata.obs.index
    adata.var_names_make_unique()
    adata_list.append(adata)

combined_adata = sc.concat(adata_list, axis=0, join="outer")
combined_adata.obs["sample_name"] = combined_adata.obs["sample_name"].astype("category")

# Randomize the rows
permuted_indices = np.random.permutation(combined_adata.n_obs)
combined_adata = combined_adata[permuted_indices, :]

# Compute QC metrics useful for downstream analysis
combined_adata.var["mito"] = combined_adata.var_names.str.upper().str.startswith(
    ("MT-")
)
combined_adata.var["ribo"] = combined_adata.var_names.str.upper().str.startswith(
    ("RPS", "RPL")
)

combined_adata.obs["mito_frac"] = np.sum(
    combined_adata[:, combined_adata.var["mito"]].X, axis=1
) / np.sum(combined_adata.X, axis=1)

sc.pp.calculate_qc_metrics(combined_adata, percent_top=[20], inplace=True)
sc.pp.calculate_qc_metrics(
    combined_adata,
    qc_vars=["ribo", "mito"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

combined_adata.write_h5ad(args.output)
