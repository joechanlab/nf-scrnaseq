process SCRAN {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/scran.sif:latest'
    publishDir "${params.outdir}/rna_scran/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path aggregation_h5ad

    output:
    val "${name}", emit: name
    path "${name}_scran.h5ad", emit: scran_h5ad

    script: // set home directory to cache basilisk
    """
    export PATH=/opt/conda/envs/scran/bin/:$PATH
    export HOME=\$PWD
    Rscript ${baseDir}/bin/scran.R \
        ${aggregation_h5ad} \
        "${name}_scran.h5ad"
    """
}
