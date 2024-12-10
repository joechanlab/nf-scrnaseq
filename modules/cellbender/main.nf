process CELLBENDER {
    label 'gpus'
    container 'us.gcr.io/broad-dsde-methods/cellbender:latest'
    containerOptions "--nv --bind ${params.mount}"
    publishDir "${params.outdir}/rna_cellbender/", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(name), path(raw_path), val(filtered_path), val(demultiplexing), val(expected_droplets)

    output:
    tuple val(name), path("${name}_cellbender.h5"), val(filtered_path), val(demultiplexing), val(expected_droplets), emit: output

    script:
    def gpu_index = task.index % params.maxForks
    if(task.executor == 'singularity')
        """
        export CUDA_VISIBLE_DEVICES=$gpu_index
        python ${baseDir}/bin/run_cellbender.py \
            ${raw_path} \
            ${name}_cellbender.h5 \
            ${expected_droplets} \
            --filtered ${filtered_path} \
            --empty_drop_training_fraction ${params.cellbender.empty_drop_training_fraction}
        """
    else
        """
        python ${baseDir}/bin/run_cellbender.py \
            ${raw_path} \
            ${name}_cellbender.h5 \
            ${expected_droplets} \
            --filtered ${filtered_path} \
            --empty_drop_training_fraction ${params.cellbender.empty_drop_training_fraction}
        """
}
