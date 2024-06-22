module load singularity/3.7.1
module load gcc/10.2.0 
module load cuda/11.7
export NXF_SINGULARITY_CACHEDIR="../work-cellbender-test/singularity"

# nextflow run ./main.nf -resume -profile singularity -params-file ./params.yml -w "../work-cellbender-test/"
nextflow run ./main.nf -resume -profile lilac -params-file ./params.yml -w "../work-cellbender-test/"