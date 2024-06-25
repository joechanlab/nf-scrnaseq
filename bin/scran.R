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
print("SCE created")
sce <- readH5AD(inpath, X_name = 'counts')

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
writeH5AD(sce, outpath, X_name = 'normcounts')
