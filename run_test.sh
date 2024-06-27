module load singularity/3.7.1
module load gcc/10.2.0
module load cuda/11.7
export NXF_SINGULARITY_CACHEDIR="/lila/data/chanjlab/wangm10/work-nf-scrnaseq/singularity/"

# nextflow run ./main.nf -resume -profile singularity -params-file ./params.yml -w "/lila/data/chanjlab/wangm10/work-nf-scrnaseq/"
# nextflow run ./main.nf -resume 0b2aab78-3380-47e5-bf3a-d2babb788b9f -profile lilac -params-file ./params.yml -w "../work-cellbender-test/"

# nextflow run ./main.nf -resume -profile singularity -params-file ../oliver_RPMN_2024/oliver_RPMN_2024_params.yml -w "/lila/data/chanjlab/wangm10/work-nf-scrnaseq/"
nextflow run ./main.nf -resume 68995d8a-02a8-478c-b01e-d1b67095de9a -profile lilac -params-file ../oliver_RPMN_2024/oliver_RPMN_2024_primary_tumor_params.yml -w "/lila/data/chanjlab/wangm10/work-nf-scrnaseq/"
# nextflow run ./main.nf -resume 5022abf3-61bc-4fe8-ab1f-18232051d5a2 -profile lilac -params-file ../oliver_RPMN_2024/oliver_RPMN_2024_params.yml -w "/lila/data/chanjlab/wangm10/work-nf-scrnaseq/"
