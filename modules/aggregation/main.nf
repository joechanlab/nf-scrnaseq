process AGGREGATION {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/doubletdetection.sif:latest'
    containerOptions "--bind ${params.mount}"
    publishDir "${params.outdir}/rna_aggregation/", mode: 'copy'
    cache 'lenient'

    input:
    val doublet_h5ads

    output:
    val "${params.experiment.name}", emit: name
    path "${params.experiment.name ? params.experiment.name + '_' : ''}aggregation*.h5ad", emit: output_h5ad

    script:
    def joined_file_paths = doublet_h5ads.join(' ')
    def agg_prefix = params.experiment.name ? params.experiment.name + '_' : ''
    def splitBy = ""
    if (params.aggregation.splitBy &&
        params.aggregation.splitBy.toString().toLowerCase() !in ['none', 'null', '']) {
        splitBy = "--split-by ${params.aggregation.splitBy}"
    }
    def metaArg = ''
    if (params.aggregation.metadata &&
        params.aggregation.metadata.toString().toLowerCase() !in ['none', 'null', '']) {
        metaArg = "--metadata ${params.aggregation.metadata}"
    }
    """
    python ${baseDir}/bin/aggregation.py \
        ${joined_file_paths} \
        --output "${agg_prefix}aggregation.h5ad" \
        ${splitBy} \
        ${metaArg}
    """
}
