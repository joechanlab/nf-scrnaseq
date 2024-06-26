# Convert 10x matrix to h5 format
# https://rdrr.io/github/MarioniLab/DropletUtils/man/write10xCounts.html
library(DropletUtils)
args <- commandArgs(trailingOnly = TRUE)
inpath = args[1]
outpath = args[2]
sce10x <- read10xCounts(inpath)
write10xCounts(outpath, sce10x@assays@data$counts, gene.id=rowData(sce10x)$ID,
    gene.symbol=rowData(sce10x)$Symbol, barcodes=colData(sce10x)$Barcode, version='3')
