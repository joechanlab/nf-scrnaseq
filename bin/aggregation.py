import argparse
import scanpy as sc
import numpy as np
from utils import get_basename_without_extension

parser = argparse.ArgumentParser(description="wrapper for concatenating the samples.")
parser.add_argument(
    "inputs", nargs="*", help="Paths to directory containing files to be concatenated."
)
parser.add_argument("output", help="The output h5ad file path.")
parser.add_argument(
    "--percent_top",
    required=False,
    type=str,
    default="50,100,200,500",
    help="Total count threshold.",
)
parser.add_argument(
    "--total_counts",
    required=False,
    type=int,
    default=500,
    help="Total count threshold.",
)
parser.add_argument(
    "--n_genes_by_counts",
    required=False,
    type=int,
    default=400,
    help="Number of genes threshold",
)
parser.add_argument(
    "--log10GenesPerUMI",
    required=False,
    type=float,
    default=0.8,
    help="log10 genes per UMI threshold",
)
parser.add_argument(
    "--mito_frac",
    required=False,
    type=float,
    default=0.2,
    help="Mitochondrial fraction threshold",
)

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

# QC
combined_adata.var["mt"] = combined_adata.var_names.str.lower().str.startswith("mt-")
combined_adata.var["ribo"] = combined_adata.var_names.str.lower().str.startswith(
    ("rps", "rpl")
)
combined_adata.var["hb"] = combined_adata.var_names.str.lower().str.contains(
    "^hb[^(p)]"
)
sc.pp.calculate_qc_metrics(
    combined_adata,
    qc_vars=["mt", "ribo", "hb"],
    inplace=True,
    log1p=True,
    percent_top=tuple(map(int, args.percent_top.split(","))),
)
combined_adata.obs["mito_frac"] = np.sum(
    combined_adata[:, combined_adata.var["mt"]].X, axis=1
) / np.sum(combined_adata.X, axis=1)
combined_adata.obs.loc[:, "log10GenesPerUMI"] = (
    np.log10(combined_adata.obs.n_genes_by_counts)
).div(np.log10(combined_adata.obs.total_counts))

sc.pp.filter_genes(combined_adata, min_cells=1)
ind1 = combined_adata.obs.total_counts > args.total_counts
ind2 = combined_adata.obs.n_genes_by_counts > args.n_genes_by_counts
ind3 = combined_adata.obs.log10GenesPerUMI > args.log10GenesPerUMI
ind4 = combined_adata.obs.mito_frac < args.mito_frac
ind = ind1.values & ind2.values & ind3.values & ind4.values
combined_adata = combined_adata[ind, :]

combined_adata.write_h5ad(args.output)
