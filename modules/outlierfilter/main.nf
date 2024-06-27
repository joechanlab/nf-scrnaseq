process OUTLIER_FILTER {
    label 'process_low'
    container 'library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest'
    publishDir "${params.outdir}/outlier_filtered/", mode: 'copy'

    input:
    path aggregation_h5ad

    output:
    path "${params.experiment.name ? params.experiment.name + '_' : ''}outlier_filtered.h5ad", emit: outlier_filtered_h5ad

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache
    python ${baseDir}/bin/outlier_filter.py \
        ${aggregation_h5ad} \
        ${params.experiment.name ? params.experiment.name + '_' : ''}outlier_filtered.h5ad
    """
}
