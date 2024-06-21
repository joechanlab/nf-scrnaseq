process CELLBENDER {
   label (params.with_gpu? 'gpus': 'process_single')
   container 'us.gcr.io/broad-dsde-methods/cellbender:latest'
   containerOptions '--nv'
   publishDir "${params.outdir}/cellbender/", mode: 'copy'

   input:
   tuple val(name), path(raw_path), path(filtered_path) // val(expected_cells), val(total_droplets_included)

   output:
   path "${name}_cellbender.h5", emit: cellbender_h5

   script:
   def gpu_index = task.index % params.maxForks
   """
   export CUDA_VISIBLE_DEVICES=$gpu_index
   python ${baseDir}/bin/run_cellbender.py \
      --raw_h5 ${raw_path} \
      --filtered_h5 ${filtered_path} \
      --output_h5 ${name}_cellbender.h5
   """
}

// cellbender remove-background \
//      --cuda \
//      --input ${path} \
//      --output ${name}_cellbender.h5 \
//      --expected-cells ${expected_cells} \
//      --total-droplets-included ${total_droplets_included}