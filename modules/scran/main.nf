process SCRAN {
   label 'process_single'
   container 'library://mamie_wang/nf-scrnaseq/scran.sif:latest'
   publishDir "${params.outdir}/scran/", mode: 'copy'

   input:
   file aggregation_h5ad
   
   output:
   path "scran.h5ad", emit: scran_h5ad

   script:
   """
   export PATH=/opt/conda/envs/scran/bin/:$PATH
   export HOME=${workDir}
   Rscript ${baseDir}/bin/scran.R \
      ${aggregation_h5ad} \
      scran.h5ad
   """
}
