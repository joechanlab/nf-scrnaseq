process CELLBENDER {
   label (params.with_gpu? 'gpus': 'process_single')
   container 'us.gcr.io/broad-dsde-methods/cellbender:latest'
   containerOptions '--nv'
   publishDir "${params.outdir}/cellbender/", mode: 'copy'

   input:
   path umi

   output:
   path "${umi.baseName}_cellbender.h5", emit: cleaned_umi

   script:
   """
   cellbender remove-background \
      --cuda \
      --input $umi \
      --output ${umi.baseName}_cellbender.h5
   """
}
