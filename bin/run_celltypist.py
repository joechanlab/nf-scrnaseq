#!/usr/bin/env python
import argparse
import os
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
original_columns = set(adata.obs.columns)
print(f"Read {args.input_h5ad}")

# Perform normalization required by celltypist
if args.normalize:
    adata.layers["scvi"] = adata.X
    adata.X = adata.layers["X_scran"]
    sc.pp.normalize_total(adata, target_sum=10**4)
    sc.pp.log1p(adata)

# Check if multiple models are specified
models_list = [model.strip() for model in args.model.split(",")]

# Initialize a list to store predictions
all_predictions = []
for model_name in models_list:
    if os.path.exists(model_name):
        print(f"Found model at {model_name}")
        model = models.Model.load(model=model_name)
        print(f"Loaded {model_name}")
    elif model_name in models.get_all_models():
        model = models.Model.load(model=model_name)
        print(f"Loaded {model_name} from celltypist's built-in models")
    else:
        print(f"Warning: Model {model_name} not found. Skipping...")
        continue
    model_adata = adata.copy()
    predictions = celltypist.annotate(
        model_adata, model=model, majority_voting=args.majority_voting
    )
    pred_adata = predictions.to_adata(
        insert_labels=True, insert_conf=True, insert_prob=True
    )
    all_predictions.append(pred_adata)

# Combine predictions if multiple models were used
if len(all_predictions) > 1:
    for idx, pred_adata in enumerate(all_predictions):
        extra_columns = set(pred_adata.obs.columns) - original_columns
        for col in extra_columns:
            new_col = f"{col}_{idx}" if col in adata.obs.columns else col
            adata.obs[new_col] = pred_adata.obs[col]
    print(f"Combined predictions from {len(all_predictions)} models")
else:
    adata = all_predictions[0].to_adata(
        insert_labels=True, insert_conf=True, insert_prob=True
    )

# Write data
adata.write_h5ad(args.output_h5ad)
print(f"Wrote {args.output_h5ad}")
