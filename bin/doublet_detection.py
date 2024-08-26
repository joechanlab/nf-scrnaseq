#!/usr/bin/env python
import os
import argparse
import scanpy as sc
import doubletdetection
from utils import anndata_from_h5, get_basename_without_extension

os.environ["NUMBA_CACHE_DIR"] = "/tmp/"

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data."
)
parser.add_argument("input_h5", help="Path to 10x h5 format.")
parser.add_argument("output_h5ad", help="The output path.")
parser.add_argument(
    "--n_iterations",
    required=False,
    default=10,
    type=int,
    help="Number of iterations to use; default is 50",
)
parser.add_argument(
    "--standard_scaling",
    required=False,
    default=True,
    type=bool,
    help="Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True",
)
parser.add_argument(
    "--p_thresh",
    required=False,
    default=1e-7,
    type=float,
    help="P-value threshold for doublet calling; default is 1e-16",
)
parser.add_argument(
    "--voter_thresh",
    required=False,
    default=0.8,
    type=float,
    help="Voter threshold for doublet calling; default is 0.5",
)
parser.add_argument(
    "--filtered_h5",
    nargs="*",
    default="",
    help="Path to cellranger filtered 10x h5 data.",
)

args = parser.parse_args()

# Read batch data
if "h5ad" in args.input_h5:
    adata_batch = sc.read_h5ad(args.input_h5)
else:
    adata_batch = anndata_from_h5(args.input_h5)

if args.filtered_h5 != "":
    print(args.filtered_h5)
    # Read sample data
    adata_sample_barcodes = []
    for idx, h5_file_path in enumerate(args.filtered_h5, start=1):
        if h5_file_path.endswith(".h5ad"):
            adata_sample = sc.read_h5ad(h5_file_path)
        else:
            adata_sample = sc.read_10x_h5(h5_file_path)
        adata_sample_barcodes += adata_sample.obs_names.to_list()

    adata_sample_barcodes = set(adata_sample_barcodes)
    adata_batch_barcodes = set(adata_batch.obs_names)

    # Find the common barcodes between sample and batch
    common_barcodes = adata_sample_barcodes.intersection(adata_batch_barcodes)

    # Filter data based on common barcodes
    adata_sample = adata_sample[adata_sample.obs_names.isin(common_barcodes)]
    adata_batch = adata_batch[adata_batch.obs_names.isin(common_barcodes)]

adata_batch = adata_batch[adata_batch.obs["cell_probability"] > 0.9]
adata_batch.var_names_make_unique()
sc.pp.filter_genes(adata_batch, min_cells=1)

# Run doublet detection
clf = doubletdetection.BoostClassifier(
    n_iters=args.n_iterations,
    clustering_algorithm="leiden",
    standard_scaling=args.standard_scaling,
    pseudocount=0.1,
    n_jobs=-1,
)

# Calculate doublets
doublets = clf.fit(adata_batch.X).predict(
    p_thresh=args.p_thresh, voter_thresh=args.voter_thresh
)
doublet_score = clf.doublet_score()

adata_batch.obs["doublet"] = doublets.astype(bool)
adata_batch.obs["doublet_score"] = doublet_score

adata_batch.obs["sample_name"] = get_basename_without_extension(args.input_h5).rsplit(
    "_", 1
)[0]
adata_batch.write_h5ad(args.output_h5ad)
