process AGGREGATION {
   label (params.with_gpu? 'gpus': 'process_single')
   container 'library://mamie_wang/nf-scrnaseq/doubletdetection.sif:latest'
   publishDir "${params.outdir}/aggregation/", mode: 'copy'

   input:
   tuple file(doublet_h5ad)
   
   output:
   path "aggregation.h5ad", emit: aggregation_h5ad

   script:
   """
   python ${baseDir}/bin/aggregation.py \
      ${doublet_h5ad} \
      --output aggregation.h5ad
   """
}
