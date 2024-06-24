
import argparse
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from utils import kneepoint
import phenograph
import scanpy as sc

parser = argparse.ArgumentParser(
    description="wrapper for postprocessing h5ad file.")
parser.add_argument('input_h5ad', help='Path to the input 10x h5ad format file.')
parser.add_argument("output_h5ad", help = "The path to the output 10x h5ad file")
parser.add_argument('--n_neighbors', type = int, default = 30, help = "Number of nearest neighbors.")
parser.add_argument('--umap_min_dist', type = float, default = 0.3, help = "Minimum distance parameter for UMAP")
parser.add_argument("--n_diffmap_components", type = int, default=2, help = "The number of diffusion map components")
parser.add_argument("--leiden_res", type = float, default = 1.8, help = "The resolution of Leiden clustering")
args = parser.parse_args()

adata = sc.read_h5ad(args.input_h5ad)

norm_df = pd.DataFrame(adata.X, index=adata.obs_names, columns = adata.var_names)

bad_genes = norm_df.columns.str.upper.str.contains(
    "^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT")
norm_df = norm_df.loc[:,~bad_genes]

'''
PCA
'''
print('Performing PCA')
n_components = 500
pca = PCA(n_components=n_components, svd_solver='randomized')
pca.fit(norm_df)

#By Kneepoint
num_components = 0
num_components = max(num_components,kneepoint(np.cumsum(pca.explained_variance_ratio_)))
print('# Components = %d' % (num_components+1))

var_explained = np.cumsum(pca.explained_variance_ratio_)[num_components]
print('Variance explained = %f' % var_explained)

pca = PCA(n_components=num_components, svd_solver='randomized')
pca_merge = pd.DataFrame(pca.fit_transform(norm_df.values),
                index=norm_df.index)
adata.obsm['X_pca'] = pca_merge.loc[adata.obs_names,:].values
adata.uns['num_components'] = num_components
adata.uns['var_explained'] = var_explained

'''
NEAREST NEIGHBORS
'''
print('Performing nearest neighbors')

sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=pca_merge.shape[1])

'''
CLUSTERING
'''
print('Phenograph Clustering')
clusters_merge, _, _ = phenograph.cluster(pca_merge, k = args.n_neighbors)
clusters_merge = pd.Series(clusters_merge, pca_merge.index)

adata.obs['phenograph'] = clusters_merge.loc[adata.obs_names].astype('str').astype('category')

'''
LEIDEN CLUSTERING
'''
print('Leiden Clustering')
sc.tl.leiden(adata, resolution = 1.8)

'''
UMAP
'''
print('Performing UMAP')
sc.tl.paga(adata, groups = 'phenograph')
sc.pl.paga(adata, plot=False)
sc.tl.umap(adata, init_pos='paga', min_dist=args.umap_min_dist)

'''
DIFFUSION MAP
'''
print('Performing Diffusion Map')
sc.tl.diffmap(adata, n_comps = args.n_diffmap_components)

'''
DEG
'''
print('Performing DEG')
sc.tl.rank_genes_groups(adata, groupby='phenograph', method='wilcoxon')

adata.write_h5ad(args.output_h5ad)
