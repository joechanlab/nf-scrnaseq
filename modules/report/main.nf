process REPORT {
    label 'process_medium'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/report/", mode: 'copy'

    input:
    path postprocessing_h5ad
    path celltypist_scvi_h5ad

    output:
    path "${params.experiment.name}_report.ipynb", emit: report_ipynb
    path "${params.experiment.name}_report.html", emit: report_html

    script:
    """
    export HOME=${workDir}
    python -m ipykernel install --user --name postprocessing
    papermill ${baseDir}/bin/QC.ipynb ${params.experiment.name}_report.ipynb -p plots ${params.report.plots}
    jupyter nbconvert --to html ${params.experiment.name}_report.ipynb
    """
}
