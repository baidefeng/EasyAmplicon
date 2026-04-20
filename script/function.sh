#!/bin/bash

# Module: function - Functional Prediction
# Function: Performs functional prediction using tools like PICRUSt, FAPROTAX, and BugBase.

# --- Default Parameters ---
DBDIR="/d/EasyMicrobiome"
THREADS=4
GROUP="Group"

# Input parameters
INPUT_FILTERED_FA_FILE="temp/filtered.fa"
INPUT_METADATA_FILE="result/metadata.txt"

# Output parameters
OUTPUT_BUGBASE_DIR="result/bugbase/"

# Intermediate files
INTERMEDIATE_GG_OTUTAB_FILE="result/gg/otutab.txt"

# --- Help Documentation ---
usage() {
    cat <<EOF
Module: function - Functional Prediction

Usage: bash easyamplicon.sh function [OPTIONS]

Options:
  -d, --dbdir DIR        Set the EasyMicrobiome database directory.
                         (Default: /d/EasyMicrobiome)
  -t, --threads INT      Set the number of threads for parallel tasks.
                         (Default: 4)
  --group STRING         Set the group column name in the metadata.
                         (Default: "Group")
  --input-filtered-fa FILE Path to the quality-filtered sequences (e.g., temp/filtered.fa).
                         (Default: ${INPUT_FILTERED_FA_FILE})
  --input-metadata-file FILE Path to the metadata file (e.g., result/metadata.txt).
                         (Default: ${INPUT_METADATA_FILE})
  --output-bugbase-dir DIR Directory for Bugbase bacterial phenotype prediction results (e.g., result/bugbase/).
                         (Default: ${OUTPUT_BUGBASE_DIR})
  -h, --help             Display this help message and exit.
EOF
    exit 0
}

# --- Parse Command Line Arguments ---
TEMP=$(getopt -o d:t:h --long dbdir:,threads:,group:,input-filtered-fa:,input-metadata-file:,output-bugbase-dir:,help -n "$0" -- "$@")
if [ $? != 0 ]; then echo "Terminating..." >&2; exit 1; fi
eval set -- "$TEMP"

while true; do
    case "$1" in
        -d|--dbdir) DBDIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --group) GROUP="$2"; shift 2 ;;
        --input-filtered-fa) INPUT_FILTERED_FA_FILE="$2"; shift 2 ;;
        --input-metadata-file) INPUT_METADATA_FILE="$2"; shift 2 ;;
        --output-bugbase-dir) OUTPUT_BUGBASE_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) echo "Internal error!"; exit 1 ;;
    esac
done

# --- Export variables for sub-scripts (if needed) ---
export db=${DBDIR}
export threads=${THREADS}
export group=${GROUP}

# --- Module Core Logic ---

## Reference-based Quantitative Feature Table

    # Align Greengenes 97% OTUs for PICRUSt/Bugbase functional prediction.
    mkdir -p "$(dirname "${INTERMEDIATE_GG_OTUTAB_FILE}")"

    # usearch alignment is faster, but if file size exceeds limit, use vsearch alignment (see Appendix 14).
    # Default: use 1 thread for <10 cores, 10 threads for >=10 cores.
    usearch -otutab "${INPUT_FILTERED_FA_FILE}" -otus ${db}/gg/97_otus.fa \
    	-otutabout "${INTERMEDIATE_GG_OTUTAB_FILE}" -threads ${threads}
    # Alignment rate 80.0%, 1 core 11m, 4 cores 3m, 10 cores 2m, memory usage 743Mb.
    head -n3 "${INTERMEDIATE_GG_OTUTAB_FILE}"

    # Statistics
    usearch -otutab_stats "${INTERMEDIATE_GG_OTUTAB_FILE}" -output "${INTERMEDIATE_GG_OTUTAB_FILE}.stat"
    cat "${INTERMEDIATE_GG_OTUTAB_FILE}.stat"

# Functional Prediction

    ### 1. Bugbase Command Line Analysis
    # Remove existing Bugbase output directory.
    # WARNING: This will delete all contents of the Bugbase output directory if it exists.

    bugbase=${db}/script/BugBase
    rm -rf "${OUTPUT_BUGBASE_DIR}"
    Rscript ${bugbase}/bin/run.bugbase.r -L ${bugbase} \
      -i "${INTERMEDIATE_GG_OTUTAB_FILE}" -m "${INPUT_METADATA_FILE}" -c ${group} -o "${OUTPUT_BUGBASE_DIR}"


      