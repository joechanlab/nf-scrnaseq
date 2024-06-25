process AGGREGATION {
    label 'process_single'
    container 'library://mamie_wang/nf-scrnaseq/doubletdetection.sif:latest'
    publishDir "${params.outdir}/aggregation/", mode: 'copy'

    input:
    tuple path(doublet_h5ad)

    output:
    path "aggregation.h5ad", emit: aggregation_h5ad

    script:
    """
    python ${baseDir}/bin/aggregation.py \
        ${doublet_h5ad} \
        aggregation.h5ad \
        --percent_top ${params.aggregation.percent_top} \
        --total_counts ${params.aggregation.total_counts} \
        --n_genes_by_counts  ${params.aggregation.n_genes_by_counts} \
        --log10GenesPerUMI ${params.aggregation.log10GenesPerUMI} \
        --mito_frac ${params.aggregation.mito_frac}
    """
}
