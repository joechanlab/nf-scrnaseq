## Introduction

A bioinformatics pipeline that preprocesses single-cell RNA-seq data. It takes a samplesheet and CellRanger raw filtered count matrix as input, performs ambient RNA correction, doublet detection, sample aggregation, batch integration.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows, where each row contains the sample name, a raw count HDF5 file and a filtered count HDF5 file from CellRanger. 

`samplesheet.csv`:
```csv
sample,raw_h5, filtered_h5
CONTROL_REP1,raw_feature_bc_matrix.h5,filtered_feature_bc_matrix.h5
```

Next, prepare a parameter YAML file that looks as follows:

`params.yml`:
```yaml
samplesheet: "./samplesheet.csv"      # path to the sample sheet
outdir: "./out/"                      # directory containing the outputs
cellbender:                           # cellbender parameters (see bin/cellbender.py)
    total_droplets_included: 50000
aggregation:                          # qc summary and filter parameters (see bin/aggregation.py)
    percent_top: "50,100,150,200" 
    total_counts: 500
    n_genes_by_counts: 400
    log10GenesPerUMI: 0.8
    mito_frac: 0.2
scvi:                                 # scvi parameters (see bin/scvi_norm.py)
    n_latent: 50
    n_top_genes: 5000
postprocessing:                       # postprocessing parameters (see bin/postprocessing.py)
   n_pca_components: 100
with_gpu: true                        # using GPU
maxForks: 2                           # maximum number of processes in parallel (e.g # of GPU)
max_memory: "6.GB"                    # memory information
max_cpus: 6                           # cpu information
```

Now, you can run the pipeline. For local run on HPC with singularity installed, execute the following command

```bash
nextflow run ./main.nf \
   -profile singularity \
   -params-file ./params.yml \
   -w ./work/
```

If you are using MSKCC lilac, you can use the pre-defined `lilac` profile that uses the LSF executor. 
```bash
nextflow run ./main.nf \
   -profile lilac \
   -params-file ./params.yml \
   -w ./work/
```
