source ~/miniconda3/etc/profile.d/conda.sh

conda activate nf-env

project_output_prefix="put_your_path_here"

nextflow run main.nf \
    --samplesheet example/samplesheet.csv \
    --region ./data/KRAS_region_hg38.bed \
    --outdir ${project_output_prefix}/nf_work

