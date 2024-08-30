process SEACELLS {
    label 'process_medium'
    conda "/usersoftware/chanj3/SEACells"
    publishDir "${params.outdir}/seacells/", mode: 'copy'

    input:
    val name
    path input_h5ad

    output:
    val "${name}", emit: name
    path "${name}_seacells.h5ad", emit: seacells_h5ad
    path "${name}_seacells_assignments", emit: seacells_assignments

    script:
    """
    export NUMBA_CACHE_DIR=\$PWD
    python ${baseDir}/bin/run_SEACells.py \
        ${input_h5ad} \
        ${name}_seacells.h5ad \
        ${name}_seacells_assignments \
        --n_SEACells ${params.seacells.n_SEACells} \
        --build_kernel_on ${params.seacells.build_kernel_on} \
        --n_waypoint_eigs ${params.seacells.n_waypoint_eigs}
    """
}
