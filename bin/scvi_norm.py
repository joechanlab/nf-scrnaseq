import os
import scvi
import scanpy as sc
import argparse

# Suppress the deprecation warning
#scvi.settings.dl_pin_memory_gpu_training = False
parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument("-i", "--input", required = True, default = "", help = "The input h5ad file")
parser.add_argument("-o", "--output", required = False, help = "The output h5ad file")
parser.add_argument('-n', '--n_latent', default=10, required=False, help='Number of latent dimensions.')

args = parser.parse_args()

adata = sc.read_h5ad(args.input)
sc.pp.log1p(adata, base=2)

adata.obs['sample_name'] = adata.obs['sample_name'].astype('category')

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=50, #5000,
    subset=True,
    flavor="seurat", 
    batch_key="sample_name"
)

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key = 'sample_name',
    categorical_covariate_keys = ['sample_name'],
    continuous_covariate_keys=["mito_frac"]
)

model = scvi.model.SCVI(adata, n_latent = args.n_latent)

model.train()

latent = model.get_latent_representation()

adata.obsm["X_scVI"] = latent

denoised = model.get_normalized_expression(adata, library_size=1e4)

adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4
)
        
adata.write_h5ad(args.output)