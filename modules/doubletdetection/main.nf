process DOUBLETDETECTION {
   label (params.with_gpu? 'gpus': 'process_single')
   container 'library://mamie_wang/nf-scrnaseq/doubletdetection.sif:latest'
   publishDir "${params.outdir}/doubletdetection/", mode: 'copy'

   input:
   tuple val(name), path(path), val(expected_cells), val(total_droplets_included)
   path cellbender_h5

   output:
   path "${name}_doubletdetection.h5ad", emit: doublet_h5ad

   script:
   """
   python ${baseDir}/bin/doublet_detection.py \
      --cellbender ${path} \
      --output ${name}_doubletdetection.h5ad
   """
}
