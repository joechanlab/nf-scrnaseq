process OUTLIER_FILTER {
    label 'process_single'
    container 'library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest'
    publishDir "${params.outdir}/outlier_filtered/", mode: 'copy'

    input:
    tuple path(aggregation_h5ad)

    output:
    path "${params.experiment.name ? params.experiment.name + '_' : ''}outlier_filtered.h5ad", emit: outlier_filtered_h5ad

    script:
    """
    python ${baseDir}/bin/outlier_filter.py \
        ${aggregation_h5ad} \
        ${params.experiment.name ? params.experiment.name + '_' : ''}outlier_filtered.h5ad \
        --percent_top ${params.outlier_filter.percent_top}
    """
}
