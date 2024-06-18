module load singularity/3.7.1
module load gcc/10.2.0 
module load cuda/11.7
export NXF_SINGULARITY_CACHEDIR="/lila/data/chanjlab/wangm10/work-cellbender-test/singularity"

nextflow run ./main.nf -resume -profile singularity -w /lila/data/chanjlab/wangm10/work-cellbender-test --outdir /lila/data/chanjlab/wangm10/results-cellbender-test --samplesheet ./data/sampleSheet2.csv