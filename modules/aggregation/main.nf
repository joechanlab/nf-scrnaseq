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
    path "${params.experiment.name ? params.experiment.name + '_' : ''}aggregation.h5ad", emit: output_h5ad

    script:
    joined_file_paths = doublet_h5ads.join(' ')
    """
    python ${baseDir}/bin/aggregation.py \
        ${joined_file_paths} \
        --output "${params.experiment.name ? params.experiment.name + '_' : ''}aggregation.h5ad"
    """
}
