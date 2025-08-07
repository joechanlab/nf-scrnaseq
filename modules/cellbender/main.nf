process CELLBENDER {
    label 'gpus'
    container 'us.gcr.io/broad-dsde-methods/cellbender:latest'
    containerOptions "--nv --bind ${params.mount}"
    publishDir "${params.outdir}/rna_cellbender/", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(name), path(raw_path), val(filtered_path), val(demultiplexing), val(expected_droplets), val(empty_drop_training_fraction)

    output:
    tuple val(name), path("${name}_cellbender.h5"), val(filtered_path), val(demultiplexing), emit: output

    script:
    def gpu_index = task.index % params.maxForks
    def has_droplets = expected_droplets && expected_droplets != "null" && expected_droplets != ""
    def has_empty_drop_training_fraction = empty_drop_training_fraction && empty_drop_training_fraction != "null" && empty_drop_training_fraction != ""

    if(task.executor == 'singularity')
        """
        export CUDA_VISIBLE_DEVICES=$gpu_index

        if [ -n "${has_droplets ? expected_droplets : ""}" ]; then
            python ${baseDir}/bin/run_cellbender.py \\
                ${raw_path} \\
                ${name}_cellbender.h5 \\
                --total_droplets_included ${expected_droplets} \\
                ${has_empty_drop_training_fraction ? "--empty_drop_training_fraction ${empty_drop_training_fraction}" : ""} \\
                --filtered ${filtered_path}
        else
            python ${baseDir}/bin/run_cellbender.py \\
                ${raw_path} \\
                ${name}_cellbender.h5 \\
                ${has_empty_drop_training_fraction ? "--empty_drop_training_fraction ${empty_drop_training_fraction}" : ""} \\
                --filtered ${filtered_path}
        fi
        """
    else
        """
        if [ -n "${has_droplets ? expected_droplets : ""}" ]; then
            python ${baseDir}/bin/run_cellbender.py \\
                ${raw_path} \\
                ${name}_cellbender.h5 \\
                --total_droplets_included ${expected_droplets} \\
                ${has_empty_drop_training_fraction ? "--empty_drop_training_fraction ${empty_drop_training_fraction}" : ""} \\
                --filtered ${filtered_path}
        else
            python ${baseDir}/bin/run_cellbender.py \\
                ${raw_path} \\
                ${name}_cellbender.h5 \\
                ${has_empty_drop_training_fraction ? "--empty_drop_training_fraction ${empty_drop_training_fraction}" : ""} \\
                --filtered ${filtered_path}
        fi
        """
}
