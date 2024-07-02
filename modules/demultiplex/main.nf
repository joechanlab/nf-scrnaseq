process DEMULTIPLEX {
    label 'process_low'
    container 'library://mamie_wang/nf-scrnaseq/doubletdetection.sif:latest'
    publishDir "${params.outdir}/demultiplex/", mode: 'copy'

    input:
    val name
    path cellbender_h5
    val filtered_path

    output:
    path "${name}_hashsolo.h5ad", emit: hashsolo_h5ad
    val ${filtered_path}, emit: filtered_path

    script:
    """
    python ${baseDir}/bin/demultiplex.py \
        ${cellbender_h5} \
        ${name}_hashsolo.h5ad
    """
}
