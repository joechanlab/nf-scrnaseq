process REPORT {
    label 'process_medium'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/report/", mode: 'copy'

    input:
    path postprocessing_h5ad
    path celltypist_scvi_h5ad

    output:
    path "${params.experiment.name}_report.html", emit: report_html

    script:
    """
    papermill ${baseDir}/bin/QC.ipynb ${params.experiment.name}_report.ipynb
    jupyter nbconvert --to html ${params.experiment.name}_report.ipynb
    """
}
