process REPORT {
    label 'process_low'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/report/", mode: 'copy'

    input:
    path postprocessing_h5ad
    path celltypist_scvi_h5ad

    output:
    path "${params.experiment.name}_QC_report.html", emit: report_html
    path "${params.experiment.name}_leiden_report.html", emit: report_html

    script:
    """
    if [ ! -d "/tmp/ipykernel" ]; then
        mkdir -p "/tmp/ipykernel"
    fi
    export HOME=/tmp/ipykernel
    python -m ipykernel install --user --name postprocessing
    papermill ${baseDir}/bin/QC.ipynb ${params.experiment.name}_QC_report.ipynb
    jupyter nbconvert --to html ${params.experiment.name}_QC_report.ipynb
    """
}
