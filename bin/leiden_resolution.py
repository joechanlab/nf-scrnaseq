import argparse
import itertools
import scanpy as sc
import pandas as pd

from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score, silhouette_score, calinski_harabasz_score
from anndata import AnnData
from typing import Optional

import contextlib

def compute_clusters_and_metrics(
    adata: AnnData,
    pca_key: str,
    r: float,
    k: int,
    random_state: int
) -> None:
    for key in [
        'leiden_q_for_r_k',
        'leiden_n_components_for_r_k',
        'leiden_n_clusters_for_r_k',
        'leiden_silhouette_coefficient_for_r_k',
        'leiden_calinski_harabasz_index_for_r_k'
    ]:
        if key not in adata.uns.keys():
            adata.uns[key] = pd.DataFrame().rename_axis(index='r', columns='k')
    clusters, graph, Q = sc.external.tl.phenograph(
        adata.obsm[pca_key],
        clustering_algo='leiden',
        k=k,
        resolution_parameter=r,
        primary_metric='correlation',
        seed=random_state,
        copy=True
    )
    adata.obs[f"leiden_r{r}_k{k}"] = pd.Categorical(clusters)
    # To allow saving to H5AD, indices cannot be floats or integers
    adata.uns['leiden_q_for_r_k'].loc[str(r), str(k)] = Q
    adata.uns['leiden_n_components_for_r_k'].loc[str(r), str(k)] = (
        connected_components(graph, directed=False, return_labels=False)
    )
    adata.uns['leiden_n_clusters_for_r_k'].loc[str(r), str(k)] = (
        adata.obs[f"leiden_r{r}_k{k}"].value_counts().count()
    )
    adata.uns['leiden_silhouette_coefficient_for_r_k'].loc[str(r), str(k)] = (
        silhouette_score(adata.obsm[pca_key], clusters)
    )
    adata.uns['leiden_calinski_harabasz_index_for_r_k'].loc[str(r), str(k)] = (
        calinski_harabasz_score(adata.obsm[pca_key], clusters)
    )

def compute_ari_matrices(
    adata: AnnData,
    rlist: list[float],
    kmin: Optional[int] = 10,
    kmax: Optional[int] = 101,
    kstep: Optional[int] = 5
) -> None:
    if 'adj_rand_index_consistency_for_r' not in adata.uns.keys():
        adata.uns['adj_rand_index_consistency_for_r'] = {}
    krange = range(kmin, kmax, kstep)
    lkrange = len(krange)
    strkrange = [str(i) for i in krange]
    for r in rlist:
        rand_score_matrix_df = pd.DataFrame(
            np.zeros((lkrange, lkrange)),
            index=strkrange,
            columns=strkrange
        )
        for i, j in itertools.combinations(strkrange, 2):
            rand_score_matrix_df.loc[i, j] = adjusted_rand_score(
                adata.obs[f"leiden_r{r}_k{i}"],
                adata.obs[f"leiden_r{r}_k{j}"]
            )
            rand_score_matrix_df.loc[j, i] = rand_score_matrix_df.loc[i, j]
            rand_score_matrix_df.loc[i, i] = 1.0
        rand_score_matrix_df.loc['100', '100'] = 1.0
        adata.uns['adj_rand_index_consistency_for_r'][str(r)] = rand_score_matrix_df

parser = argparse.ArgumentParser(description="wrapper for postprocessing h5ad file.")
parser.add_argument("input", help="Path to the input 10x h5ad format file.")
parser.add_argument("output", help="The path to the output directory")
parser.add_argument(
    "--res_range", type=str, default="0.25,0.5,0.75,1,1.25", help="Resolution to check."
)
args = parser.parse_args()

if not args.use_scvi:
    dim_red = 'X_pca'
else:
    dim_red = "X_scvi"

adata = sc.read_h5ad(args.input)
rlist = args.res_range
for r in rlist:
    for k in range(10, 101, 5):
        with contextlib.redirect_stdout(None):
            compute_clusters_and_metrics(adata, 'X_pca', r, k, random_state=123456)

compute_ari_matrices(adata, rlist=rlist)

print("Write H5AD")
adata.write_h5ad(args.output_h5ad)

print(f"Wrote H5AD as {args.output_h5ad}")
