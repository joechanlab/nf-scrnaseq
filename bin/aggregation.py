import os
import re
import argparse
import scanpy as sc
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="wrapper for concatenating the samples.")
parser.add_argument(
    "inputs", nargs="*", help="Paths to directory containing files to be concatenated."
)
parser.add_argument("--output", required=True, help="The output h5ad file path.")
parser.add_argument(
    "--split-by", help="Column in .obs to split the data.", default=None
)
parser.add_argument(
    "--metadata", required=False, default="None", help="Metadata to be added to obs."
)
args = parser.parse_args()

adata_list = []

for idx, h5_file_path in enumerate(args.inputs, start=1):
    print(h5_file_path)
    adata = sc.read(h5_file_path)
    adata.obs.index = adata.obs["sample_name"].astype(str) + "_" + adata.obs.index
    adata.var_names_make_unique()
    adata_list.append(adata)

combined_adata = sc.concat(adata_list, axis=0, join="outer")
combined_adata.obs["sample_name"] = combined_adata.obs["sample_name"].astype("category")

# Randomize the rows
permuted_indices = np.random.permutation(combined_adata.n_obs)
combined_adata = combined_adata[permuted_indices, :]

if args.metadata and args.metadata != "None":
    print("Reading the metadata...")
    metadata = pd.read_csv(args.metadata)
    new_cols = [x for x in metadata.columns if x not in combined_adata.obs.columns]
    intersect_cols = [x for x in metadata.columns if x in combined_adata.obs.columns]
    metadata = combined_adata.obs.merge(metadata, how="left", on=intersect_cols)
    assert metadata.shape[0] == combined_adata.shape[0]
    metadata.index = combined_adata.obs.index
    for new_col in new_cols:
        s = metadata[new_col]
        non_na = s.dropna()
        unique_types = {type(v) for v in non_na}
        if (
            s.dtype == bool
            or unique_types <= {bool, np.bool_, int, np.int64}
            and set(non_na.unique()) <= {0, 1, True, False}
        ):
            combined_adata.obs[new_col] = s.astype("boolean")
        elif s.dtype == "object":
            combined_adata.obs[new_col] = s.astype(str)
        else:
            combined_adata.obs[new_col] = s

if not args.split_by:
    print(f"Writing combined aggregation to: {args.output}")
    combined_adata.write_h5ad(args.output)
else:
    missing_col = args.split_by not in combined_adata.obs.columns
    if missing_col:
        raise KeyError(
            f"The following split-by columns are missing in .obs: {args.split_by}. "
            f"Available columns include: {list(combined_adata.obs.columns)}"
        )
    # remove nan columns
    combined_adata.obs[args.split_by] = combined_adata.obs[args.split_by].replace(
        ["nan", "NaN", "None", ""], np.nan
    )
    combined_adata = combined_adata[
        combined_adata.obs[args.split_by].notnull(), :
    ].copy()
    base, ext = os.path.splitext(args.output)
    if ext == "":
        ext = ".h5ad"

    def _sanitize(val: str) -> str:
        """Make a safe string for filenames."""
        return re.sub(r"[^A-Za-z0-9]+", "-", str(val)).strip("-")

    # Group by requested obs columns and save per-group
    grouped = combined_adata.obs.groupby(args.split_by)
    print(f"Splitting aggregation by: {args.split_by}")
    for group_vals, indices in grouped.groups.items():
        if not isinstance(group_vals, tuple):
            group_vals = (group_vals,)

        suffix_parts = []
        for val in group_vals:
            suffix_parts.append(f"{_sanitize(val)}")
        suffix = "__" + "__".join(suffix_parts)

        out_path = f"{base}{suffix}{ext}"
        print(f"Writing group {group_vals} to: {out_path}")
        combined_adata[indices, :].write_h5ad(out_path)
