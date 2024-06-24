#!/usr/bin/env python
import os
import argparse
import scanpy as sc
import doubletdetection
from utils import anndata_from_h5

parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
parser.add_argument('input_h5', default='', required=True, help='Path to 10x h5 format.')
parser.add_argument("output", help = "The output directory; default is current working directory")
parser.add_argument("--n_iterations", required = False, default = 10, type = int, help = "Number of iterations to use; default is 50")
parser.add_argument("--standard_scaling", required = False, default = True, type = bool, help = "Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True")
parser.add_argument("--p_thresh", required = False, default = 1e-7, type = float, help = "P-value threshold for doublet calling; default is 1e-16")
parser.add_argument("--voter_thresh", required = False, default = 0.8, type = float, help = "Voter threshold for doublet calling; default is 0.5")
parser.add_argument('--filtered_h5', default='', nargs = "*", help='Path to cellranger filtered 10x h5 data.')

args = parser.parse_args()

# Read batch data
adata_batch = anndata_from_h5(args.input_h5)

if args.filtered_h5 != '':
    # Read sample data
    adata_sample_barcodes = []
    for idx, h5_file_path in enumerate(args.inputs, start=1):
        adata_sample = sc.read_10x_h5(args.filtered_h5)
        # Extract sample barcodes
        adata_sample_barcodes.append(adata_sample.obs_names)

    adata_sample_barcodes = set(adata_sample_barcodes)
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