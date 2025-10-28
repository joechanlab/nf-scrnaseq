export NXF_SINGULARITY_CACHEDIR="/usersoftware/chanj3/singularity/"

nextflow run ./main.nf -profile iris -params-file ./params.yml -w "/scratch/chanj3/wangm10/work"
