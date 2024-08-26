import argparse
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import statsmodels.formula.api as smf

from anndata import AnnData
from typing import Optional, Callable

# Arguments
parser = argparse.ArgumentParser(
    description="Outlier detection from transcriptomic data."
)
parser.add_argument("input_h5ad", help="Path to 10x h5ad format.")
parser.add_argument("output_h5ad", help="The output path.")
parser.add_argument(
    "--sample_col", default="sample_name", help="The sample column name."
)
args = parser.parse_args()


# Definitions of helper functions
def calculate_featcount_dist(
    adata: AnnData, *, key_added: Optional[str] = "featcount_dist"
) -> None:
    """
    Calculates feature-count distances and stores them in adata.obs[key_added].
    Feature-count distance is the difference between the observed and expected
    log no. of features given log total counts.
    From Germain et al. 2020 (DOI: 10.1186/s13059-020-02136-7).
    """
    mod = smf.ols("log1p_n_genes_by_counts ~ log1p_total_counts", adata.obs).fit()
    pred = mod.predict()
    adata.obs[key_added] = adata.obs["log1p_n_genes_by_counts"] - pred


def calculate_group_featcount_dist(
    adata: AnnData, *, group_key: str, key_added: Optional[str] = None
) -> None:
    """
    Calculates group-wise feature-count distances for groups in
    adata.obs[group_key] and stores them in adata.obs[key_added]
    """
    if key_added is None:
        key_added = f"featcount_dist_{group_key}"
    groups = adata.obs[group_key].unique().tolist()
    for group in groups:
        tmp_adata = adata[adata.obs[group_key] == group].copy()
        calculate_featcount_dist(tmp_adata, key_added="tmp_fcd")
        adata.obs.loc[adata.obs[group_key] == group, key_added] = tmp_adata.obs[
            "tmp_fcd"
        ]


def mad_outlier(
    adata: AnnData, *, metric: str, nmads_upper: float, nmads_lower: float
) -> pd.Series:
    """
    Marks cell as outlier (i.e. returns True) if it is more than
    nmads_upper/nmads_lower median absolute deviations above/below
    the median of the metric.
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads_lower * sp.stats.median_abs_deviation(M)) | (
        M > np.median(M) + nmads_upper * sp.stats.median_abs_deviation(M)
    )
    return outlier


def default_filter(adata: AnnData, group_key: str) -> np.ndarray:
    """
    Default filter adapted from Germain et al. 2020 (DOI: 10.1186/s13059-020-02136-7)
    and Heumos et al. 2023 (DOI: 10.1038/s41576-023-00586-w). Outliers must meet at
    least two of the five given criteria:
    1. log1p_total_counts < median - 3 MADs OR log1p_total_counts > median + 5 MADs
    2. log1p_n_genes_by_counts < median - 3 MADs OR log1p_n_genes_by_counts > median + 5 MADs
    3. pct_counts_in_top_20_genes < median - 5 MADs OR pct_counts_in_top_20_genes > median + 5 MADs
    4. featcount_dist_<group> < median - 5 MADs OR featcount_dist_<group> > median + 5 MADs
    5. pct_counts_mito > median + 2.5 MADs AND pct_counts_mito > 8%
    """
    condition = (
        np.sum(
            [
                mad_outlier(
                    adata, metric="log1p_total_counts", nmads_lower=3, nmads_upper=5
                ),
                mad_outlier(
                    adata,
                    metric="log1p_n_genes_by_counts",
                    nmads_lower=3,
                    nmads_upper=5,
                ),
                mad_outlier(
                    adata,
                    metric="pct_counts_in_top_20_genes",
                    nmads_lower=5,
                    nmads_upper=5,
                ),
                mad_outlier(
                    adata,
                    metric=f"featcount_dist_{group_key}",
                    nmads_lower=5,
                    nmads_upper=5,
                ),
                (
                    mad_outlier(
                        adata,
                        metric="pct_counts_mito",
                        nmads_lower=np.inf,
                        nmads_upper=2.5,
                    )
                    & (adata.obs["pct_counts_mito"] > 8)
                ),
            ],
            axis=0,
        )
        >= 2
    )
    return condition


def designate_outliers(
    adata: AnnData,
    *,
    condition: Callable[[AnnData, str], np.ndarray],
    group_key: str,
    key_added: Optional[str] = "outlier",
) -> None:
    """
    Marks group-wise outliers according to condition given as a function
    of adata and group_key, and saves them to adata.obs[key_added]
    """
    groups = adata.obs[group_key].unique().tolist()
    for group in groups:
        adata.obs.loc[adata.obs[group_key] == group, key_added] = condition(
            adata[adata.obs[group_key] == group], group_key
        )
    adata.obs[key_added] = adata.obs[key_added].astype(bool)


# Load data
adata = sc.read_h5ad(args.input_h5ad)

# Compute QC metrics useful for downstream analysis
adata.var["mito"] = adata.var_names.str.upper().str.startswith(("MT-"))
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))

adata.obs["mito_frac"] = np.sum(adata[:, adata.var["mito"]].X, axis=1) / np.sum(
    adata.X, axis=1
)

sc.pp.calculate_qc_metrics(adata, percent_top=[20], inplace=True)
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["ribo", "mito"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

calculate_group_featcount_dist(adata, group_key=args.sample_col)

# Mark outliers
designate_outliers(adata, condition=default_filter, group_key=args.sample_col)

# Filter genes
sc.pp.filter_genes(adata, min_cells=1)

# Save results
adata.write_h5ad(args.output_h5ad)
