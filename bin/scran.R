# Load libraries
suppressMessages({library(scran)
                  library(zellkonverter)
                  library(scater)
                 })
args <- commandArgs(trailingOnly = TRUE)

# Load DGC matrix file and metadata
inpath <- args[1]
outpath <- args[2]

print("Data loaded")

# Create single cell experiment
sce <- readH5AD(inpath)

print("SCE created")
assays(sce)[['counts']] <- assays(sce)[['X']]
assays(sce)[['X']] <- NULL

# Calculate stats
qcstats <- perCellQCMetrics(sce)
qcfilter <- quickPerCellQC(qcstats)
sce <- sce[,!qcfilter$discard]

print("QC filter")
summary(qcfilter$discard)

# Compute clusters
clusters <- quickCluster(sce)

# Compute factors
sce <- computeSumFactors(sce, clusters=clusters)

print("Size factors")
summary(sizeFactors(sce))

# Normalize
sce <- logNormCounts(sce, log=F)

# save the file
assays(sce)[['X']] <- assays(sce)[['normcounts']]
assays(sce)[['counts']] <- NULL
assays(sce)[['normcounts']] <- NULL
writeH5AD(sce, outpath)
