process SCRAN {
    label 'process_medium'
    container 'library://mamie_wang/nf-scrnaseq/scran.sif:latest'
    publishDir "${params.outdir}/scran/", mode: 'copy'

    input:
    file aggregation_h5ad

    output:
    path "${params.experiment.name ? params.experiment.name + '_' : ''}scran.h5ad", emit: scran_h5ad

    script: // set home directory to cache basilisk
    """
    export PATH=/opt/conda/envs/scran/bin/:$PATH
    export HOME=${workDir}
    Rscript ${baseDir}/bin/scran.R \
        ${aggregation_h5ad} \
        "${params.experiment.name ? params.experiment.name + '_' : ''}scran.h5ad"
    """
}
