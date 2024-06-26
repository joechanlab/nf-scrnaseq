process DOUBLETDETECTION {
    label 'process_single'
    container 'library://mamie_wang/nf-scrnaseq/doubletdetection.sif:latest'
    publishDir "${params.outdir}/doubletdetection/", mode: 'copy'

    input:
    val name
    path cellbender_h5
    val filtered_path

    output:
    path "${name}_doubletdetection.h5ad", emit: doublet_h5ad

    script:
    """
    python ${baseDir}/bin/doublet_detection.py \
        ${cellbender_h5} \
        ${name}_doubletdetection.h5ad \
        --filtered_h5 ${filtered_path}
    """
}
