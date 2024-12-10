process POSTPROCESSING {
    label 'process_medium'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/rna_postprocessing/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path scvi_h5ad

    output:
    val "${name}", emit: name
    path "${name}_postprocessing.h5ad", emit: postprocessing_h5ad
    path "${name}_postprocessing_scvi.h5ad", emit: postprocessing_scvi_h5ad

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/postprocessing.py \
        ${scvi_h5ad} \
        ${name}_postprocessing.h5ad \
        --n_pca_components ${params.postprocessing.n_pca_components} \
        --metadata ${params.postprocessing.metadata} \
        --leiden_res ${params.postprocessing.leiden_res}
    python ${baseDir}/bin/postprocessing.py \
        ${scvi_h5ad} \
        ${name}_postprocessing_scvi.h5ad \
        --use_scvi \
        --metadata ${params.postprocessing.metadata} \
        --leiden_res ${params.postprocessing.leiden_res}
    """
}
