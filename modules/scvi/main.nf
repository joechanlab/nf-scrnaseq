process SCVI {
    label 'gpus'
    container "docker://scverse/scvi-tools:py3.11-cu12-tutorials-latest"
    containerOptions '--nv'
    publishDir "${params.outdir}/scvi/", mode: 'copy'

    input:
    path scran_h5ad

    output:
    path "${params.experiment.name ? params.experiment.name + '_' : ''}scvi.h5ad", emit: scvi_h5ad

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/scvi_norm.py \
        ${scran_h5ad} \
        ${params.experiment.name ? params.experiment.name + '_' : ''}scvi.h5ad \
        --n_latent ${params.scvi.n_latent} \
        --n_top_genes ${params.scvi.n_top_genes}
    """
}
