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
    def has_droplets = expected_droplets && expected_droplets != "null" && expected_droplets != ""

    if(task.executor == 'singularity')
        """
        export CUDA_VISIBLE_DEVICES=$gpu_index
        if [ -n "${has_droplets ? expected_droplets : ""}" ]; then
            python ${baseDir}/bin/run_cellbender.py \\
                ${raw_path} \\
                ${name}_cellbender.h5 \\
                ${expected_droplets} \\
                --filtered ${filtered_path}
        else
            python ${baseDir}/bin/run_cellbender.py \\
                ${raw_path} \\
                ${name}_cellbender.h5 \\
                --filtered ${filtered_path}
        fi
        """
    else
        """
        if [ -n "${has_droplets ? expected_droplets : ""}" ]; then
            python ${baseDir}/bin/run_cellbender.py \\
                ${raw_path} \\
                ${name}_cellbender.h5 \\
                ${expected_droplets} \\
                --filtered ${filtered_path}
        else
            python ${baseDir}/bin/run_cellbender.py \\
                ${raw_path} \\
                ${name}_cellbender.h5 \\
                --filtered ${filtered_path}
        fi
        """
}
