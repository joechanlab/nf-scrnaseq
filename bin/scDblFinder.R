suppressMessages({library(zellkonverter)
    library(BiocParallel)
    library(scDblFinder)
    })
args <- commandArgs(trailingOnly = TRUE)

INPUT_PATH <- args[1]
OUTPUT_PATH <- args[2]
N_CORES <- args[3]
SEED <- 18591124

# Load data; X should be untransformed counts
sce <- readH5AD(INPUT_PATH, X_name = "counts")

# Set seeds for replicability
bp <- MulticoreParam(workers = N_CORES, RNGseed = SEED)
set.seed(SEED)

# Run scDblFinder
sce <- scDblFinder(sce, BPPARAM = bp)

writeH5AD(sce, OUTPUT_PATH)
