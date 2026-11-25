#!/usr/bin/env bash
# Set up the software for bamrc_mutation_curation:
#   1. a conda environment with Nextflow and Python
#   2. the bam-readcount and SnpEff Singularity images
#
# Usage:
#   bash setup.sh [-n ENV_NAME] [-d IMAGE_DIR]
#
#   -n  conda environment name           (default: nf-env)
#   -d  directory to save the .sif files (default: ./images)
#
# After it finishes, put the printed image paths into your config
# (bamrc_sif, snpeff_sif). See README.md.

set -euo pipefail

ENV_NAME="nf-env"
IMAGE_DIR="./images"

BAMRC_IMAGE="docker://mgibio/bam-readcount:1.0.0"
SNPEFF_IMAGE="docker://lihao5657/snpeff:4.3t"

while getopts "n:d:h" opt; do
    case $opt in
        n) ENV_NAME="$OPTARG" ;;
        d) IMAGE_DIR="$OPTARG" ;;
        h) sed -n '2,14p' "$0"; exit 0 ;;
        *) sed -n '2,14p' "$0"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------- checks
if command -v mamba >/dev/null 2>&1; then
    CONDA=mamba
elif command -v conda >/dev/null 2>&1; then
    CONDA=conda
else
    echo "ERROR: conda (or mamba) not found on PATH." >&2
    exit 1
fi

if command -v singularity >/dev/null 2>&1; then
    SINGULARITY=singularity
elif command -v apptainer >/dev/null 2>&1; then
    SINGULARITY=apptainer
else
    echo "ERROR: singularity/apptainer not found on PATH (try 'module load singularity')." >&2
    exit 1
fi

# ---------------------------------------------------------------- conda env
echo "==> Conda environment: ${ENV_NAME}"
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "    already exists, skipping (remove it with 'conda env remove -n ${ENV_NAME}' to rebuild)"
else
    "${CONDA}" create -y -n "${ENV_NAME}" \
        -c conda-forge -c bioconda \
        "nextflow>=23.10" "python>=3.8"
fi

# ---------------------------------------------------------------- images
mkdir -p "${IMAGE_DIR}"
IMAGE_DIR="$(cd "${IMAGE_DIR}" && pwd)"

pull() {
    local out="$1" src="$2"
    if [ -s "${out}" ]; then
        echo "    ${out} already exists, skipping"
    else
        "${SINGULARITY}" pull "${out}" "${src}"
    fi
}

echo "==> Pulling images into ${IMAGE_DIR}"
pull "${IMAGE_DIR}/bam-readcount_1.0.0.sif" "${BAMRC_IMAGE}"
pull "${IMAGE_DIR}/snpeff_4.3t.sif"         "${SNPEFF_IMAGE}"

# ---------------------------------------------------------------- next steps
cat <<EOF

Setup complete.

Set these in your config (e.g. site.config, passed with -c):

    params {
        bamrc_sif  = '${IMAGE_DIR}/bam-readcount_1.0.0.sif'
        snpeff_sif = '${IMAGE_DIR}/snpeff_4.3t.sif'
    }

Then activate the environment before running the pipeline:

    conda activate ${ENV_NAME}
EOF
