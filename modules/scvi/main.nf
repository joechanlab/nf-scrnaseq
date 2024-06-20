process SCVI {
   label (params.with_gpu? 'gpus': 'process_single')
   container "library://mamie_wang/nf-scrnaseq/scvi.sif:latest"
   publishDir "${params.outdir}/scvi/", mode: 'copy'

   input:
   file scran_h5ad
   
   output:
   path "scvi.h5ad", emit: scvi_h5ad

   script:
   """
   export NUMBA_CACHE_DIR=/tmp/numba_cache
   python ${baseDir}/bin/scvi_norm.py \
      --input ${scran_h5ad} \
      --output scvi.h5ad
   """
}
