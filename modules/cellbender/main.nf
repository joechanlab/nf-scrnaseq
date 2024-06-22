process CELLBENDER {
   label 'gpus'
   container 'us.gcr.io/broad-dsde-methods/cellbender:latest'
   containerOptions '--nv'
   publishDir "${params.outdir}/cellbender/", mode: 'copy'

   input:
   tuple val(name), path(raw_path), path(filtered_path) 

   output:
   path "${name}_cellbender.h5", emit: cellbender_h5

   script:
   """
   python ${baseDir}/bin/run_cellbender.py \
      --raw_h5 ${raw_path} \
      --filtered_h5 ${filtered_path} \
      --output_h5 ${name}_cellbender.h5 \
      --total_droplets_included ${params.cellbender.total_droplets_included}
   """
}

// def gpu_index = task.index % params.maxForks
// export CUDA_VISIBLE_DEVICES=$gpu_index