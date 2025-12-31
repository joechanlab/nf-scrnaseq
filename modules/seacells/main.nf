process SEACELLS {
    label 'process_medium'
    conda "/usersoftware/chanj3/SEACells"
    publishDir "${params.outdir}/rna_seacells/", mode: 'copy'
    cache 'lenient'

    input:
    val name
    path input_h5ad

    output:
    val "${name}", emit: name
    path "*_seacells.h5ad", emit: seacells_h5ad
    path "*_seacells_assignments", emit: seacells_assignments

    script:
    def full  = input_h5ad.baseName
    def base = full.contains("__") ? full.split("__", 2)[1] : full
    """
    export NUMBA_CACHE_DIR=\$PWD
    /usersoftware/chanj3/SEACells/bin/python ${baseDir}/bin/run_SEACells.py \
        ${input_h5ad} \
        ${base}_seacells.h5ad \
        ${base}_seacells_assignments \
        --n_SEACells ${params.seacells.n_SEACells} \
        --build_kernel_on ${params.seacells.build_kernel_on} \
        --n_waypoint_eigs ${params.seacells.n_waypoint_eigs}
    """
}
