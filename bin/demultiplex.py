import scanpy.external as sce
from utils import anndata_from_h5
import argparse

parser = argparse.ArgumentParser(description="Demultiplexing.")
parser.add_argument("input", help="Paths to the raw h5ad file.")
parser.add_argument("--output", required=True, help="The output h5ad file path.")

args = parser.parse_args()

adata = anndata_from_h5(args.input)

multiplexed = adata.var["feature_type"] == "Multiplexing Capture"

if sum(multiplexed) > 0:
    htos = adata[:, multiplexed].var_names.tolist()

    # Copy HTO data to obs
    for hto in htos:
        adata.obs[hto] = adata[:, adata.var_names == hto].X.toarray()

    # Use flat prior for pre-QC data
    sce.pp.hashsolo(
        adata,
        cell_hashing_columns=htos,
        priors=[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        inplace=True,
    )

    # Remove HTO data from obs
    for hto in htos:
        del adata.obs[hto]

# Save file with all HashSolo info
adata.write_h5ad(args.output)

# Filter cells and features
# adata = adata[~adata.obs['Classification'].isin(['Doublet', 'Negative']), :].copy()
# adata = adata[:, adata.var['feature_type'] != 'Multiplexing Capture'].copy()

# # Delete HashSolo info and save cleaned data
# hashsolo_obs = [
#     'most_likely_hypothesis', 'cluster_feature',
#     'negative_hypothesis_probability', 'singlet_hypothesis_probability',
#     'doublet_hypothesis_probability'
# ]

# for key in hashsolo_obs:
#     del adata.obs[key]

# adata.write_h5ad(FILTERED_OUTPUT_PATH)
