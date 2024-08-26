process SCVI {
    label 'gpus'
    container "docker://scverse/scvi-tools:py3.11-cu12-tutorials-latest"
    containerOptions '--nv'
    publishDir "${params.outdir}/scvi/", mode: 'copy'

    input:
    val name
    path scran_h5ad

    output:
    val "${name}", emit: name
    path "${name}_scvi.h5ad", emit: scvi_h5ad

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/scvi_norm.py \
        ${scran_h5ad} \
        ${name}_scvi.h5ad \
        --n_latent ${params.scvi.n_latent} \
        --n_top_genes ${params.scvi.n_top_genes}
    """
}
