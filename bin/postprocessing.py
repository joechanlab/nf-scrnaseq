import argparse
import numpy as np
import pandas as pd
import numpy.matlib
from sklearn.decomposition import PCA
import scanpy as sc


def kneepoint(vec):
    curve = [1 - x for x in vec]
    nPoints = len(curve)
    allCoord = np.vstack((range(nPoints), curve)).T
    np.array([range(nPoints), curve])
    firstPoint = allCoord[0]
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec**2))
    vecFromFirst = allCoord - firstPoint
    scalarProduct = np.sum(
        vecFromFirst * numpy.matlib.repmat(lineVecNorm, nPoints, 1), axis=1
    )
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = np.sqrt(np.sum(vecToLine**2, axis=1))
    idxOfBestPoint = np.argmax(distToLine)
    return idxOfBestPoint


def RunPCA(cts, var_threshold, n_components):
    pca = PCA(n_components=n_components, svd_solver="randomized")
    pca.fit(cts)
    num_components = 0
    num_components = max(
        num_components, kneepoint(np.cumsum(pca.explained_variance_ratio_))
    )
    num_components = max(
        num_components,
        np.where(np.cumsum(pca.explained_variance_ratio_) > var_threshold)[0][0],
    )
    var_explained = np.cumsum(pca.explained_variance_ratio_)[num_components]
    print("# Components = %d" % (num_components + 1))
    print("Variance explained = %f" % var_explained)
    return pca, num_components, var_explained


parser = argparse.ArgumentParser(description="wrapper for postprocessing h5ad file.")
parser.add_argument("input_h5ad", help="Path to the input 10x h5ad format file.")
parser.add_argument("output_h5ad", help="The path to the output 10x h5ad file")
parser.add_argument(
    "--n_neighbors", type=int, default=30, help="Number of nearest neighbors."
)
parser.add_argument(
    "--umap_min_dist",
    type=float,
    default=0.3,
    help="Minimum distance parameter for UMAP",
)
parser.add_argument(
    "--n_pca_components", type=int, default=500, help="Number of PCA components"
)
parser.add_argument(
    "--n_diffmap_components",
    type=int,
    default=3,
    help="The number of diffusion map components",
)
parser.add_argument(
    "--leiden_res", type=float, default=1.8, help="The resolution of Leiden clustering"
)
args = parser.parse_args()

adata = sc.read_h5ad(args.input_h5ad)

norm_df = pd.DataFrame(
    adata.X.toarray(), index=adata.obs_names, columns=adata.var_names
)
bad_genes = norm_df.columns.str.upper().str.contains(
    "^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT"
)
norm_df = norm_df.loc[:, ~bad_genes]

"""
PCA
"""
print("Performing PCA")
pca = PCA(n_components=args.n_pca_components, svd_solver="randomized")
pca.fit(norm_df)

# By Kneepoint
num_components = 0
num_components = max(
    num_components, kneepoint(np.cumsum(pca.explained_variance_ratio_))
)
print("# Components = %d" % (num_components + 1))

var_explained = np.cumsum(pca.explained_variance_ratio_)[num_components]
print("Variance explained = %f" % var_explained)

pca = PCA(n_components=num_components, svd_solver="randomized")
pca_merge = pd.DataFrame(pca.fit_transform(norm_df.values), index=norm_df.index)
adata.obsm["X_pca"] = pca_merge.loc[adata.obs_names, :].values
adata.uns["num_components"] = num_components
adata.uns["var_explained"] = var_explained

"""
NEAREST NEIGHBORS
"""
print("Performing nearest neighbors")
sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=pca_merge.shape[1])

"""
LEIDEN CLUSTERING
"""
print("Leiden Clustering")
sc.tl.leiden(adata, resolution=args.leiden_res)

"""
UMAP
"""
print("Performing UMAP")
sc.tl.paga(adata, groups="leiden")
sc.pl.paga(adata, plot=False)
sc.tl.umap(adata, init_pos="paga", min_dist=args.umap_min_dist)

"""
DIFFUSION MAP
"""
print("Performing Diffusion Map")
sc.tl.diffmap(adata, n_comps=args.n_diffmap_components)

"""
DEG
"""
print("Performing DEG")
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

"""
Save H5AD
"""
print("Write H5AD")
adata.write_h5ad(args.output_h5ad)

print(f"Wrote H5AD as {args.output_h5ad}")
