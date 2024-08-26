import argparse
import scanpy as sc
import numpy as np

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
    adata.obs.index = adata.obs["sample_name"] + "_" + adata.obs.index
    adata.var_names_make_unique()
    adata_list.append(adata)

combined_adata = sc.concat(adata_list, axis=0, join="outer")
combined_adata.obs["sample_name"] = combined_adata.obs["sample_name"].astype("category")

# Randomize the rows
permuted_indices = np.random.permutation(combined_adata.n_obs)
combined_adata = combined_adata[permuted_indices, :]

combined_adata.write_h5ad(args.output)
