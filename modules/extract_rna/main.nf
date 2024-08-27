process EXTRACT_RNA {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/muon.sif:latest'
    containerOptions "--bind ${params.mount}"
    publishDir "${params.outdir}/extract_rna/", mode: 'copy'

    input:
    tuple val(name), path(raw_path), path(filtered_path), val(demultiplexing), val(expected_droplets)

    output:
    tuple val(name), path("${name}_raw_RNA.h5ad"), path("${name}_filtered_RNA.h5ad"), val(demultiplexing), val(expected_droplets), emit: output

    script:
    """
    export PATH=/opt/conda/envs/muon/bin/:$PATH
    python ${baseDir}/bin/extract_rna.py \
        ${raw_path} \
        ${name}_raw_RNA.h5ad
    python ${baseDir}/bin/extract_rna.py \
        ${filtered_path} \
        ${name}_filtered_RNA.h5ad
    """
}
