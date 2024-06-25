import scvi
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data."
)
parser.add_argument(
    "-i", "--input", required=True, default="", help="The input h5ad file"
)
parser.add_argument("-o", "--output", required=False, help="The output h5ad file")
parser.add_argument(
    "--n_latent",
    default=10,
    type=int,
    required=False,
    help="Number of latent dimensions.",
)
parser.add_argument(
    "--n_top_genes", default=100, type=int, required=False, help="Number of top genes."
)

args = parser.parse_args()

adata = sc.read_h5ad(args.input)
adata.layers["without_log"] = adata.X
adata.raw = adata

sc.pp.log1p(adata, base=2)

adata.obs["sample_name"] = adata.obs["sample_name"].astype("category")

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=args.n_top_genes,
    subset=True,
    flavor="seurat",
    batch_key="sample_name",
)

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="sample_name",
    categorical_covariate_keys=["sample_name"],
    continuous_covariate_keys=["mito_frac"],
)

model = scvi.model.SCVI(adata, n_latent=args.n_latent)

model.train()

latent = model.get_latent_representation()

adata.obsm["X_scVI"] = latent

denoised = model.get_normalized_expression(adata, library_size=1e4)

adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=10e4)

adata.write_h5ad(args.output)
