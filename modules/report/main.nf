process REPORT {
    label 'process_medium'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/rna_report/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path postprocessing_h5ad
    path celltypist_scvi_h5ad

    output:
    path "${name}_report.ipynb", emit: report_ipynb
    path "${name}_report.html", emit: report_html

    script:
    def plotsParam = params.report.plot
    def plotsArg = ''
    if (!plotsParam || plotsParam.toString().trim().toLowerCase() in ['none','null','']) {
        plotsArg = "-r plots None"
    } else {
        plotsArg = "-p plots ${plotsParam}"
    }
    """
    export HOME=\$PWD
    python -m ipykernel install --user --name postprocessing
    papermill ${baseDir}/bin/QC.ipynb ${name}_report.ipynb ${plotsArg}
    jupyter nbconvert --to html ${name}_report.ipynb
    """
}
