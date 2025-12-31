import scvi
import scanpy as sc
import argparse
import torch

S_GENES = [
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM2",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDCA7",
    "DTL",
    "PRIM1",
    "UHRF1",
    "MLF1IP",
    "HELLS",
    "RFC2",
    "RPA2",
    "NASP",
    "RAD51AP1",
    "GMNN",
    "WDR76",
    "SLBP",
    "CCNE2",
    "UBR7",
    "POLD3",
    "MSH2",
    "ATAD2",
    "RAD51",
    "RRM2",
    "CDC45",
    "CDC6",
    "EXO1",
    "TIPIN",
    "DSCC1",
    "BLM",
    "CASP8AP2",
    "USP1",
    "CLSPN",
    "POLA1",
    "CHAF1B",
    "BRIP1",
    "E2F8",
]

G2M_GENES = [
    "HMGB2",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "NUF2",
    "CKS1B",
    "MKI67",
    "TMPO",
    "CENPF",
    "TACC3",
    "FAM64A",
    "SMC4",
    "CCNB2",
    "CKAP2L",
    "CKAP2",
    "AURKB",
    "BUB1",
    "KIF11",
    "ANP32E",
    "TUBB4B",
    "GTSE1",
    "KIF20B",
    "HJURP",
    "CDCA3",
    "HN1",
    "CDC20",
    "TTK",
    "CDC25C",
    "KIF2C",
    "RANGAP1",
    "NCAPD2",
    "DLGAP5",
    "CDCA2",
    "CDCA8",
    "ECT2",
    "KIF23",
    "HMMR",
    "AURKA",
    "PSRC1",
    "ANLN",
    "LBR",
    "CKAP5",
    "CENPE",
    "CTCF",
    "NEK2",
    "G2E3",
    "GAS2L3",
    "CBX5",
    "CENPA",
]

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
    "--sample_key",
    default="sample_name",
    help="Sample key to be used for highly variable genes.",
)
parser.add_argument(
    "--categorical_covariates",
    nargs="*",
    default=None,
    help="Space-separated list of obs columns to use as categorical covariates "
    "(e.g. CYW transduction).",
)
parser.add_argument(
    "--continuous_covariates",
    nargs="*",
    default=["mito_frac", "S_score", "G2M_score"],
    help="Space-separated list of obs columns to use as continuous covariates "
    "(default: mito_frac).",
)

args = parser.parse_args()

torch.set_float32_matmul_precision("high")

adata = sc.read_h5ad(args.input)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

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
s_genes = [g for g in S_GENES if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
cat_covs = []
if args.categorical_covariates is not None:
    cat_covs = [c for c in args.categorical_covariates if c in adata_hvg.obs.columns]

cont_covs = []
if args.continuous_covariates is not None:
    cont_covs = [c for c in args.continuous_covariates if c in adata_hvg.obs.columns]

scvi.model.SCVI.setup_anndata(
    adata_hvg,
    batch_key="sample_name",
    layer="counts",
    categorical_covariate_keys=cat_covs if len(cat_covs) > 0 else None,
    continuous_covariate_keys=cont_covs if len(cont_covs) > 0 else None,
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
