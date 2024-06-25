import os
import scvi
import scanpy as sc
import argparse
import torch

parser = argparse.ArgumentParser(
    description="wrapper for running SCVI on transcriptomic data.")
parser.add_argument("input", help="The input h5ad file")
parser.add_argument("output", help="The output h5ad file")
parser.add_argument("--n_latent", default=10, type=int, required=False, help="Number of latent dimensions.")
parser.add_argument("--n_top_genes", default=1000, type=int, required=False, help="Number of top genes.")

args = parser.parse_args()

torch.set_float32_matmul_precision("high")

adata = sc.read_h5ad(args.input)
adata.layers["X_scran"] = adata.X
adata.raw = adata
sc.pp.log1p(adata, base=2)

adata.obs["sample_name"] = adata.obs["sample_name"].astype("category")

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=args.n_top_genes,
    flavor="seurat_v3", 
    batch_key="sample_name"
)

adata_hvg = adata[:, adata.var["highly_variable"]].copy()

scvi.model.SCVI.setup_anndata(
    adata_hvg,
    batch_key="sample_name",
    continuous_covariate_keys=["mito_frac"]
)

model = scvi.model.SCVI(adata_hvg, n_latent=args.n_latent)

scvi.settings.seed = 18591124
model.train(
    early_stopping=True,
    early_stopping_monitor='reconstruction_loss_validation'
)

scvi.settings.seed = 18591124
latent = model.get_latent_representation()

scvi.settings.seed = 18591124
denoised = model.get_normalized_expression(library_size=1e4)

scvi.settings.seed = 18591124
integrated = model.get_normalized_expression(
    transform_batch=adata.obs["sample_name"].unique().tolist(),
    library_size=1e4
)

adata.obsm["X_scVI"] = latent
adata.obsm["scvi_normalized"] = denoised
adata.obsm["scvi_normalized_integrated"] = integrated
        
adata.write_h5ad(args.output)
