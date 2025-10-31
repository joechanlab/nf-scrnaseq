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
    // Choose a GPU (customize to your scheme)
    def gpu_index = task.index % (params.maxForks ?: 1)

    def _toStr = { it?.toString()?.trim() }
    def expStr   = _toStr(expected_droplets)
    def emptyStr = _toStr(empty_drop_training_fraction)

    def expFlag   = (expStr && expStr != 'null')   ? "--total_droplets_included ${expStr}" : ''
    def emptyFlag = (emptyStr && emptyStr != 'null') ? "--empty_drop_training_fraction ${emptyStr}" : ''

    // For Singularity we often need to export CUDA_VISIBLE_DEVICES explicitly
    def maybeCudaExport = (task.executor == 'singularity') ? "export CUDA_VISIBLE_DEVICES=${gpu_index}\n" : ""

    """
    ${maybeCudaExport}python ${baseDir}/bin/run_cellbender.py \\
        ${raw_path} \\
        ${name}_cellbender.h5 \\
        ${expFlag} \\
        ${emptyFlag} \\
        --filtered ${filtered_path}
    """
}
