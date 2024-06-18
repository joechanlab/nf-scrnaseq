#!/usr/bin/env python
import os
import argparse
import tables
import numpy as np
import scipy.sparse as sp
import anndata
import scanpy as sc
import doubletdetection
from typing import Dict, Optional

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument('-c', '--cellbender', default='', required=True, help='Path to 10x h5 format.')
parser.add_argument("-o", "--output", required = True, default = os.getcwd(), help = "The output directory; default is current working directory")
parser.add_argument('-i', '--umi', default='', required=False, help='Path to 10x h5 format.')
parser.add_argument("-n", "--n_iterations", required = False, default = 10, type = int, help = "Number of iterations to use; default is 50")
parser.add_argument("-s", "--standard_scaling", required = False, default = True, type = bool, help = "Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True")
parser.add_argument("-t", "--p_thresh", required = False, default = 1e-7, type = float, help = "P-value threshold for doublet calling; default is 1e-16")
parser.add_argument("-v", "--voter_thresh", required = False, default = 0.8, type = float, help = "Voter threshold for doublet calling; default is 0.5")

args = parser.parse_args()

def get_basename_without_extension(path):
    basename = os.path.basename(path) 
    basename_without_extension = os.path.splitext(basename)[0]
    return basename_without_extension


def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary.

    Args:
        file: The h5 file

    Returns:
        Dictionary containing all the information from the h5 file
    """
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(file: str,
                    analyzed_barcodes_only: bool = True) -> anndata.AnnData:
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count matrix.
            True to load a limited set of barcodes: only those analyzed by the
            algorithm. This allows relevant latent variables to be loaded
            properly into adata.obs and adata.obsm, rather than adata.uns.

    Returns:
        anndata.AnnData: The anndata object, populated with inferred latent variables
            and metadata.

    """

    d = dict_from_h5(file)
    X = sp.csc_matrix((d.pop('data'), d.pop('indices'), d.pop('indptr')),
                      shape=d.pop('shape')).transpose().tocsr()

    # check and see if we have barcode index annotations, and if the file is filtered
    barcode_key = [k for k in d.keys() if (('barcode' in k) and ('ind' in k))]
    if len(barcode_key) > 0:
        max_barcode_ind = d[barcode_key[0]].max()
        filtered_file = (max_barcode_ind >= X.shape[0])
    else:
        filtered_file = True

    if analyzed_barcodes_only:
        if filtered_file:
            # filtered file being read, so we don't need to subset
            print('Assuming we are loading a "filtered" file that contains only cells.')
            pass
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        elif 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        else:
            print('Warning: analyzed_barcodes_only=True, but the key '
                  '"barcodes_analyzed_inds" or "barcode_indices_for_latents" '
                  'is missing from the h5 file. '
                  'Will output all barcodes, and proceed as if '
                  'analyzed_barcodes_only=False')

    # Construct the anndata object.
    adata = anndata.AnnData(X=X,
                            obs={'barcode': d.pop('barcodes').astype(str)},
                            var={'gene_name': (d.pop('gene_names') if 'gene_names' in d.keys()
                                               else d.pop('name')).astype(str)},
                            dtype=X.dtype)
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_name', inplace=True)

    # For CellRanger v2 legacy format, "gene_ids" was called "genes"... rename this
    if 'genes' in d.keys():
        d['id'] = d.pop('genes')

    # For purely aesthetic purposes, rename "id" to "gene_id"
    if 'id' in d.keys():
        d['gene_id'] = d.pop('id')

    # If genomes are empty, try to guess them based on gene_id
    if 'genome' in d.keys():
        if np.array([s.decode() == '' for s in d['genome']]).all():
            if '_' in d['gene_id'][0].decode():
                print('Genome field blank, so attempting to guess genomes based on gene_id prefixes')
                d['genome'] = np.array([s.decode().split('_')[0] for s in d['gene_id']], dtype=str)

    # Add other information to the anndata object in the appropriate slot.
    _fill_adata_slots_automatically(adata, d)

    # Add a special additional field to .var if it exists.
    if 'features_analyzed_inds' in adata.uns.keys():
        adata.var['cellbender_analyzed'] = [True if (i in adata.uns['features_analyzed_inds'])
                                            else False for i in range(adata.shape[1])]
    elif 'features_analyzed_inds' in adata.var.keys():
        adata.var['cellbender_analyzed'] = [True if (i in adata.var['features_analyzed_inds'].values)
                                            else False for i in range(adata.shape[1])]

    if analyzed_barcodes_only:
        for col in adata.obs.columns[adata.obs.columns.str.startswith('barcodes_analyzed')
                                     | adata.obs.columns.str.startswith('barcode_indices')]:
            try:
                del adata.obs[col]
            except Exception:
                pass
    else:
        # Add a special additional field to .obs if all barcodes are included.
        if 'barcodes_analyzed_inds' in adata.uns.keys():
            adata.obs['cellbender_analyzed'] = [True if (i in adata.uns['barcodes_analyzed_inds'])
                                                else False for i in range(adata.shape[0])]
        elif 'barcodes_analyzed_inds' in adata.obs.keys():
            adata.obs['cellbender_analyzed'] = [True if (i in adata.obs['barcodes_analyzed_inds'].values)
                                                else False for i in range(adata.shape[0])]

    return adata

def _fill_adata_slots_automatically(adata, d):
    """Add other information to the adata object in the appropriate slot."""

    for key, value in d.items():
        try:
            if value is None:
                continue
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == adata.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == adata.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print('Unable to load data into AnnData: ', key, value, type(value))

# Read batch data
adata_batch = anndata_from_h5(args.cellbender)

if args.umi != '':
    # Read sample data
    adata_sample = sc.read_10x_mtx(args.umi)
    sample_name = get_basename_without_extension(args.umi)

    # Extract sample barcodes
    adata_sample_barcodes = set(adata_sample.obs_names)

    # Extract batch barcodes
    adata_batch_barcodes = set(adata_batch.obs_names)

    # Find the common barcodes between sample and batch
    common_barcodes = adata_sample_barcodes.intersection(adata_batch_barcodes)

    # Filter data based on common barcodes
    adata_sample = adata_sample[adata_sample.obs_names.isin(common_barcodes)]
    adata_batch = adata_batch[adata_batch.obs_names.isin(common_barcodes)]

adata_batch = adata_batch[adata_batch.obs['cell_probability'] > 0.9]
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
doublets = clf.fit(adata_batch.X).predict(p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)
doublet_score = clf.doublet_score()

adata_batch.obs["doublet"] = doublets
adata_batch.obs["doublet_score"] = doublet_score

adata_batch.write_h5ad(args.output)

# Save convergence plot with sample name
# f = doubletdetection.plot.convergence(clf, save = f"{args.outdir}/{sample_name}.convergence_test.pdf", show=True, p_thresh=1e-7, voter_thresh=0.8)

# sc.pp.normalize_total(adata_batch)
# sc.pp.log1p(adata_batch)
# sc.pp.highly_variable_genes(adata_batch)
# sc.tl.pca(adata_batch)
# sc.pp.neighbors(adata_batch)
# sc.tl.umap(adata_batch)

# Save UMAP plot with sample name
# umap_plot_file = f"{args.outdir}/{sample_name}_umap_doublets.png"
# sc.pl.umap(adata_batch, color=["doublet", "doublet_score"], save = umap_plot_file)

# Save violin plot with sample name
# violin_plot_file = f"{args.outdir}/{sample_name}_violin_doublets.png"
# sc.pl.violin(adata_batch, "doublet_score", save = violin_plot_file)

# Save threshold plot with sample name
# threshold_plot_file = f"{args.outdir}/{sample_name}_threshold_test.pdf"
# f3 = doubletdetection.plot.threshold(clf, save = threshold_plot_file, show=True, p_step=6)