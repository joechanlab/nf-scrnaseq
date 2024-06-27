process CELLBENDER {
    label 'gpus'
    container 'us.gcr.io/broad-dsde-methods/cellbender:latest'
    containerOptions '--nv'
    publishDir "${params.outdir}/cellbender/", mode: 'copy'

    input:
    tuple val(name), path(raw_path), val(filtered_path)

    output:
    val "${name}", emit: name
    path "${name}_cellbender.h5", emit: cellbender_h5
    val "${filtered_path}", emit: filtered_path

    script:
    def gpu_index = task.index % params.maxForks
    if(task.executor == 'lsf')
        """
        python ${baseDir}/bin/run_cellbender.py \
            ${raw_path} \
            ${name}_cellbender.h5 \
            ${params.cellbender.total_droplets_included} \
            --filtered ${filtered_path}
        """
    else
        """
        export CUDA_VISIBLE_DEVICES=$gpu_index
        python ${baseDir}/bin/run_cellbender.py \
            ${raw_path} \
            ${name}_cellbender.h5 \
            ${params.cellbender.total_droplets_included} \
            --filtered ${filtered_path}
        """
}
