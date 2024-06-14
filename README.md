## Introduction

A bioinformatics pipeline that preprocesses single-cell RNA-seq data. It takes a samplesheet and CellRanger raw filtered count matrix as input, performs ambient RNA correction, doublet detection, sample aggregation, batch integration.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,raw_h5
CONTROL_REP1,raw_feature_bc_matrix.h5
```

Each row represents a raw count matrix file 

Now, you can run the pipeline using:

```bash
nextflow run mamie/nf-scrnaseq \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```