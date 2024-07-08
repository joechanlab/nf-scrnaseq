#!/usr/bin/env python
import argparse
import scanpy as sc
import celltypist
from celltypist import models

parser = argparse.ArgumentParser(
    description="Run celltypist for cell type annotation. "
)
parser.add_argument("input_h5ad", help="Path to h5ad file.")
parser.add_argument("output_h5ad", help="The output path to h5ad file.")
parser.add_argument(
    "--model",
    default="Human_Lung_Atlas.pkl",
    type=str,
    help="Type of model; default is Human_Lung_Atlas.pkl",
)
parser.add_argument(
    "--majority_voting",
    type=bool,
    default=False,
    action=argparse.BooleanOptionalAction,
    help="Whether to use majority voting classifier.",
)
parser.add_argument(
    "--use_gpu",
    type=bool,
    default=False,
    action=argparse.BooleanOptionalAction,
    help="Whether to use GPU for annotaiton.",
)
parser.add_argument(
    "--normalize",
    type=bool,
    default=False,
    action=argparse.BooleanOptionalAction,
    help="Whether to normalize the data.",
)

args = parser.parse_args()

# Read data
adata = sc.read_h5ad(args.input_h5ad)
print(f"Read {args.input_h5ad}")

# Perform normalization required by celltypist
if args.normalize:
    adata.layers["scvi"] = adata.X
    adata.X = adata.layers["X_scran"]
    sc.pp.normalize_total(adata, target_sum=10**4)
    sc.pp.log1p(adata)

# Download and download the model
models.download_models(model=args.model)
print(f"Downloaded {args.model} at {models.models_path}")

# Load the model
model = models.Model.load(model=args.model)
print(f"Loaded {args.model}")

# prediction
predictions = celltypist.annotate(
    adata, model=model, majority_voting=args.majority_voting
)
print("Made predictions")

# Write data
adata = predictions.to_adata(insert_labels=True, insert_conf=True, insert_prob=True)
adata.write_h5ad(args.output_h5ad)
print(f"Wrote {args.output_h5ad}")
