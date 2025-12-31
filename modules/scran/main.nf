process SCRAN {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/scran.sif:latest'
    publishDir "${params.outdir}/rna_scran/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path outlier_filter_h5ad

    output:
    val "${name}", emit: name
    path "*_scran.h5ad", emit: scran_h5ad

    script: // set home directory to cache basilisk
    def full = outlier_filter_h5ad.baseName
    def base = (full.contains("__")) ? full.split("__", 2)[1] : full
    """
    export PATH=/opt/conda/envs/scran/bin/:$PATH
    export HOME=\$PWD
    Rscript ${baseDir}/bin/scran.R \
        ${outlier_filter_h5ad} \
        "${base}_scran.h5ad"
    """
}
