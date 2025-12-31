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
    path "*_report.ipynb", emit: report_ipynb
    path "*_report.html", emit: report_html

    script:
    def full  = postprocessing_h5ad.baseName
    def base = full.contains("__") ? full.split("__", 2)[1] : full
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
    papermill ${baseDir}/bin/QC.ipynb ${base}_report.ipynb ${plotsArg}
    jupyter nbconvert --to html ${base}_report.ipynb
    """
}
