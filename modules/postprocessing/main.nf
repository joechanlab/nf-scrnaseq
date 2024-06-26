process POSTPROCESSING {
    label 'process_single'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/postprocessing/", mode: 'copy'

    input:
    path scvi_h5ad

    output:
    path "postprocessing.h5ad", emit: postprocessing_h5ad

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache
    python ${baseDir}/bin/postprocessing.py \
        ${scvi_h5ad} \
        postprocessing.h5ad \
        --n_pca_components ${params.postprocessing.n_pca_components} || true
    """
}
