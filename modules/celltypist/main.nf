process CELLTYPIST {
    label 'gpus'
    container 'quay.io/teichlab/celltypist:latest'
    containerOptions '--nv'
    publishDir "${params.outdir}/celltypist/", mode: 'copy'

    input:
    path postprocessing_scvi_h5ad

    output:
    path "${params.experiment.name}_celltypist_scvi.h5ad", emit: celltypist_scvi_h5ad

    script:
    def gpu_index = task.index % params.maxForks
    if(task.executor == 'singularity')
        """
        export CUDA_VISIBLE_DEVICES=\$PWD
        export NUMBA_CACHE_DIR=\$PWD
        export MPLCONFIGDIR=\$PWD
        python ${baseDir}/bin/run_celltypist.py \
            ${postprocessing_scvi_h5ad} \
            ${params.experiment.name}_celltypist_scvi.h5ad \
            --model ${params.celltypist.model} \
            --majority_voting \
            --use_gpu \
            --normalize
        """
    else
        """
        export NUMBA_CACHE_DIR=\$PWD
        export MPLCONFIGDIR=\$PWD
        python ${baseDir}/bin/run_celltypist.py \
            ${postprocessing_scvi_h5ad} \
            ${params.experiment.name}_celltypist_scvi.h5ad \
            --model ${params.celltypist.model} \
            --majority_voting \
            --use_gpu \
            --normalize
        """
}
