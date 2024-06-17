process CELLBENDER {
   label (params.with_gpu? 'gpus': 'process_single')
   container 'us.gcr.io/broad-dsde-methods/cellbender:latest'
   containerOptions '--nv'
   publishDir "${params.outdir}/cellbender/", mode: 'copy'

   input:
   tuple val(name), path(path), val(expected_cells)

   output:
   path "${name}_cellbender.h5", emit: cellbender_h5

   script:
   """
   cellbender remove-background \
      --cuda \
      --input ${path} \
      --expected-cells ${expected_cells} \
      --output ${name}_cellbender.h5
   """
}