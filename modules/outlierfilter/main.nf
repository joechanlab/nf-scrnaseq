process OUTLIER_FILTER {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest'
    publishDir "${params.outdir}/rna_outlier_filtered/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path aggregation_h5ad

    output:
    val "${name}", emit: name
    path "${name}_outlier_filtered.h5ad", emit: outlier_filtered_h5ad

    script:
    if (params.remove_outliers)
        """
        export NUMBA_CACHE_DIR=\$PWD
        python ${baseDir}/bin/outlier_filter.py \
            ${aggregation_h5ad} \
            ${name}_outlier_filtered.h5ad \
            --remove_outliers
        """
    else
        """
        export NUMBA_CACHE_DIR=\$PWD
        python ${baseDir}/bin/outlier_filter.py \
            ${aggregation_h5ad} \
            ${name}_outlier_filtered.h5ad
        """
}
