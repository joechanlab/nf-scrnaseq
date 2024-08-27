process DOUBLETDETECTION {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/doubletdetection.sif:latest'
    containerOptions "--bind ${params.mount}"
    publishDir "${params.outdir}/doubletdetection/", mode: 'copy'

    input:
    val name
    path cellbender_h5
    path filtered_path
    val demultiplexing

    output:
    val "${name}", emit: name
    path "${name}_doubletdetection.h5ad", emit: doublet_h5ad

    script:
    if(demultiplexing == 'true')
        """
        export NUMBA_CACHE_DIR=\$PWD
        python ${baseDir}/bin/demultiplex.py \
            ${cellbender_h5} \
            --output ${name}_hashsolo.h5ad
        python ${baseDir}/bin/doublet_detection.py \
            ${name}_hashsolo.h5ad \
            ${name}_doubletdetection.h5ad \
            --filtered_h5 ${filtered_path}
        """
    else if (params.remove_doublets)
        """
        export NUMBA_CACHE_DIR=\$PWD
        python ${baseDir}/bin/doublet_detection.py \
            ${cellbender_h5} \
            ${name}_doubletdetection.h5ad \
            --filtered_h5 ${filtered_path} \
            --remove_doublets
        """
    else
        """
        export NUMBA_CACHE_DIR=\$PWD
        python ${baseDir}/bin/doublet_detection.py \
            ${cellbender_h5} \
            ${name}_doubletdetection.h5ad \
            --filtered_h5 ${filtered_path}
        """
}
