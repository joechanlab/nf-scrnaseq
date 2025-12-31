process SCVI {
    label 'gpus'
    container "docker://scverse/scvi-tools:py3.11-cu12-tutorials-latest"
    containerOptions "--nv --bind ${params.mount} --no-home"
    publishDir "${params.outdir}/rna_scvi/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path scran_h5ad

    output:
    val "${name}", emit: name
    path "*_scvi.h5ad", emit: scvi_h5ad

    script:
    def full  = scran_h5ad.baseName
    def base = full.contains("__") ? full.split("__", 2)[1] : full
    def sampleKeyArg = ''
    if (params.scvi.sample_key &&
        params.scvi.sample_key.toString().toLowerCase() !in ['none', 'null', '']) {
        sampleKeyArg = "--sample_key ${params.scvi.sample_key}"
    }

    def catCovArg = ''
    if (params.scvi.categorical_covariates) {
        // allow either list or comma-separated string in nextflow.config
        def cats = params.scvi.categorical_covariates
        if (cats instanceof String) {
            cats = cats.split(/[,\s]+/).findAll { it }    // split on comma/space, drop empties
        }
        if (cats && cats.size() > 0) {
            catCovArg = "--categorical_covariates " + cats.join(' ')
        }
    }

    def contCovArg = ''
    if (params.scvi.continuous_covariates) {
        def conts = params.scvi.continuous_covariates
        if (conts instanceof String) {
            conts = conts.split(/[,\s]+/).findAll { it }
        }
        if (conts && conts.size() > 0) {
            contCovArg = "--continuous_covariates " + conts.join(' ')
        }
    }
    """
    export PYTHONNOUSERSITE=1
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/scvi_norm.py \
        ${scran_h5ad} \
        ${base}_scvi.h5ad \
        --n_latent ${params.scvi.n_latent} \
        --n_top_genes ${params.scvi.n_top_genes} \
        ${sampleKeyArg} \
        ${catCovArg} \
        ${contCovArg}
    """
}
