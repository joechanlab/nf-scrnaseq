process EXTRACT_RNA {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/muon.sif:latest'
    containerOptions "--bind ${params.mount}"
    publishDir "${params.outdir}/extract_rna/", mode: 'copy'
    cache 'lenient'

    input:
    tuple val(name), path(raw_path), val(filtered_path), val(demultiplexing), val(expected_droplets)

    output:
    tuple val(name), path("${name}_raw_RNA.h5ad"), path("${name}_filtered_RNA.h5ad"), val(demultiplexing), val(expected_droplets), emit: output
    path "${name}_raw_ATAC.h5ad", emit: atac_raw
    path "${name}_filtered_ATAC.h5ad", emit: atac_filtered

    script:
    """
    export PATH=/opt/conda/envs/muon/bin/:$PATH
    python ${baseDir}/bin/extract_rna.py ${raw_path} ${name}_raw_RNA.h5ad ${name}_raw_ATAC.h5ad
    python ${baseDir}/bin/extract_rna.py ${filtered_path} ${name}_filtered_RNA.h5ad ${name}_filtered_ATAC.h5ad
    """
}
