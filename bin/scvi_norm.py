import scvi
import scanpy as sc
import argparse
import torch
import pandas as pd

parser = argparse.ArgumentParser(
    description="wrapper for running SCVI on transcriptomic data."
)
parser.add_argument("input", help="The input h5ad file")
parser.add_argument("output", help="The output h5ad file")
parser.add_argument(
    "--n_latent",
    default=10,
    type=int,
    required=False,
    help="Number of latent dimensions.",
)
parser.add_argument(
    "--n_top_genes", default=1000, type=int, required=False, help="Number of top genes."
)
parser.add_argument(
    "--metadata", required=False, default="None", help="Metadata to be added to obs."
)
parser.add_argument(
    "--sample_key",
    default="sample_name",
    help="Sample key to be used for highly variable genes.",
)

args = parser.parse_args()

torch.set_float32_matmul_precision("high")

adata = sc.read_h5ad(args.input)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

if not args.metadata == "None":
    print("Reading the metadata...")
    metadata = pd.read_csv(args.metadata)
    new_cols = [x for x in metadata.columns if x not in adata.obs.columns]
    intersect_cols = [x for x in metadata.columns if x in adata.obs.columns]
    metadata = adata.obs.merge(metadata, how="left", on=intersect_cols)
    assert metadata.shape[0] == adata.shape[0]
    metadata.index = adata.obs.index
    for new_col in new_cols:
        adata.obs[new_col] = metadata[new_col]

adata.layers["X_scran"] = adata.X.copy()
sc.pp.log1p(adata, base=2)
adata = adata[adata.obs[args.sample_key].notnull(), :].copy()
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=args.n_top_genes,
    layer="counts",
    flavor="seurat_v3",
    batch_key=args.sample_key,
    span=0.5,
)

adata_hvg = adata[:, adata.var["highly_variable"]].copy()

scvi.model.SCVI.setup_anndata(
    adata_hvg,
    batch_key="sample_name",
    layer="counts",
    continuous_covariate_keys=["mito_frac"],
)

model = scvi.model.SCVI(adata_hvg, n_latent=args.n_latent)

scvi.settings.seed = 18591124
model.train(
    early_stopping=True, early_stopping_monitor="reconstruction_loss_validation"
)

scvi.settings.seed = 18591124
latent = model.get_latent_representation()

scvi.settings.seed = 18591124
denoised = model.get_normalized_expression(library_size=1e4, return_numpy=True)

scvi.settings.seed = 18591124
integrated = model.get_normalized_expression(
    transform_batch=adata.obs["sample_name"].unique().tolist(),
    library_size=1e4,
    return_numpy=True,
)

adata.obsm["X_scVI"] = latent
adata.obsm["scvi_normalized"] = denoised
adata.obsm["scvi_normalized_integrated"] = integrated

adata.write_h5ad(args.output)
