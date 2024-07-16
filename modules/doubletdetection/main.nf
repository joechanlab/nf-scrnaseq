process DOUBLETDETECTION {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/doubletdetection.sif:latest'
    containerOptions '--bind /lila:/lila'
    publishDir "${params.outdir}/doubletdetection/", mode: 'copy'

    input:
    val name
    path cellbender_h5
    val filtered_path
    val demultiplexing

    output:
    path "${name}_doubletdetection.h5ad", emit: doublet_h5ad

    script:
    if(demultiplexing)
        """
        export NUMBA_CACHE_DIR=/tmp/numba_cache
        python ${baseDir}/bin/demultiplex.py \
            ${cellbender_h5} \
            --output ${name}_hashsolo.h5ad
        python ${baseDir}/bin/doublet_detection.py \
            ${name}_hashsolo.h5ad \
            ${name}_doubletdetection.h5ad \
            --filtered_h5 ${filtered_path}
        """
    else
        """
        export NUMBA_CACHE_DIR=/tmp/numba_cache
        python ${baseDir}/bin/doublet_detection.py \
            ${cellbender_h5} \
            ${name}_doubletdetection.h5ad \
            --filtered_h5 ${filtered_path}
        """
}
