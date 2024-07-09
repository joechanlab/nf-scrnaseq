process LEIDEN_RESOLUTION {
    label 'process_low'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/leiden_resolution/", mode: 'copy'

    input:
    path postprocessing_scvi_h5ad

    output:
    path "${params.experiment.name ? params.experiment.name + '_' : ''}leiden_resolution.h5ad", emit: leiden_resolution_h5ad

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache
    python ${baseDir}/bin/leiden_resolution.py \
        ${postprocessing_scvi_h5ad} \
        ${params.experiment.name ? params.experiment.name + '_' : ''}leiden_resolution.h5ad \
        --use_scvi
    papermill ${baseDir}/bin/leiden_resolution_.ipynb ${params.experiment.name}_leiden_report.ipynb
    jupyter nbconvert --to html ${params.experiment.name}_leiden_report.ipynb
    """
}
