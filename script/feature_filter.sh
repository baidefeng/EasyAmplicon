#!/bin/bash

# Module: feature_filter - Feature table normalization (assumes taxonomy/nonBac already done)

# Define general functions
check_env_var() {
    for var in "$@"; do
        if [ -z "${!var}" ]; then
            echo "Error: Environment variable '${var}' is not set." >&2
            exit 1
        fi
    done
}

# --- Default Parameters ---
# EasyMicrobiome root (contains script/, etc.). Override per machine:
#   export EASY_MICROBIOME=/path/to/EasyMicrobiome
#   or: bash easyamplicon.sh feature_filter --dbdir /path/to/EasyMicrobiome ...
DBDIR="${EASY_MICROBIOME:-/d/EasyMicrobiome}"
OUTPUT_DIR="result"

# Input parameters
INPUT_OTUTAB_FILE="${OUTPUT_DIR}/otutab.txt"
INPUT_OTUS_FILE="${OUTPUT_DIR}/otus.fa"

# Output parameters
OUTPUT_OTUTAB_STAT_FILE="${OUTPUT_DIR}/otutab.stat"
OUTPUT_ALPHA_DIVERSITY_FILE="${OUTPUT_DIR}/alpha/vegan.txt"
OUTPUT_RAREFIED_OTUTAB_FILE="${OUTPUT_DIR}/otutab_rare.txt"

# Intermediate temporary files (normalization only)
TEMP_INPUT_OTUTAB_FILE="${OUTPUT_DIR}/temp_input_otutab.txt"
DEPTH=10000 # A reasonable default, or user must provide

# --- Help Documentation ---
usage() {
    cat <<EOF
Module: feature_filter - Feature Table Filtering and Normalization

Usage: bash easyamplicon.sh feature_filter [OPTIONS]

Options:
  -d, --dbdir DIR        Set the EasyMicrobiome database directory.
                         (Default: ${DBDIR})
  --output-dir DIR       Main results output directory.
                         (Default: ${OUTPUT_DIR})
  --input-otutab-file FILE Input feature table (OTU/ASV table).
                         (Default: ${INPUT_OTUTAB_FILE})
  --input-otus-file FILE Input feature sequences.
                         (Default: ${INPUT_OTUS_FILE})
  --output-otutab-stat-file FILE Output file for feature table statistics.
                         (Default: ${OUTPUT_OTUTAB_STAT_FILE})
  --output-alpha-diversity-file FILE Output file for Alpha diversity indices.
                         (Default: ${OUTPUT_ALPHA_DIVERSITY_FILE})
  --output-rarefied-otutab-file FILE Output file for rarefied feature table.
                         (Default: ${OUTPUT_RAREFIED_OTUTAB_FILE})
  --depth INT            Set the sampling depth for rarefaction.
                         (Default: ${DEPTH})
  -h, --help             Display this help message and exit.
EOF
    exit 0
}

# --- Parse Command Line Arguments ---
TEMP=$(getopt -o d:h --long dbdir:,output-dir:,input-otutab-file:,input-otus-file:,output-otutab-stat-file:,output-alpha-diversity-file:,output-rarefied-otutab-file:,depth:,help -n "$0" -- "$@")
if [ $? != 0 ]; then echo "Terminating..." >&2; exit 1; fi
eval set -- "$TEMP"

while true; do
    case "$1" in
        -d|--dbdir) DBDIR="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --input-otutab-file) INPUT_OTUTAB_FILE="$2"; shift 2 ;;
        --input-otus-file) INPUT_OTUS_FILE="$2"; shift 2 ;;
        --output-otutab-stat-file) OUTPUT_OTUTAB_STAT_FILE="$2"; shift 2 ;;
        --output-alpha-diversity-file) OUTPUT_ALPHA_DIVERSITY_FILE="$2"; shift 2 ;;
        --output-rarefied-otutab-file) OUTPUT_RAREFIED_OTUTAB_FILE="$2"; shift 2 ;;
        --depth) DEPTH="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) echo "Internal error!"; exit 1 ;;
    esac
done

# --- Export variables for sub-scripts (if needed) ---
export db=${DBDIR}

# Check necessary environment variables
check_env_var "db" "OUTPUT_DIR" "INPUT_OTUTAB_FILE" "INPUT_OTUS_FILE" "OUTPUT_OTUTAB_STAT_FILE" "OUTPUT_ALPHA_DIVERSITY_FILE" "OUTPUT_RAREFIED_OTUTAB_FILE" "DEPTH"

# Helpers
require_file() {
    local f="$1"
    local what="${2:-file}"
    if [ ! -f "$f" ]; then
        echo "Error: missing ${what}: $f" >&2
        echo "  Download link (TODO): <PUT_DOWNLOAD_LINK_HERE>" >&2
        exit 1
    fi
}

# We assume taxonomy/nonBac filtering has already been done upstream (feature_table/feature_pacbio/feature_nanopore).
if [ ! -s "${INPUT_OTUTAB_FILE}" ]; then
    echo "Error: input feature table is missing/empty: ${INPUT_OTUTAB_FILE}" >&2
    exit 1
fi
require_file "${db}/script/otutab_rare.R" "R script otutab_rare.R"

### Equal-depth Sampling Normalization

    echo "--- Equal-depth Sampling Normalization ---"
    # Normalize by subsample


    if [[ -z "${DEPTH}" ]]; then
        echo "Error: Sampling depth is not set. Please specify using --depth parameter, or ensure it can be automatically determined." >&2
        exit 1
    fi

    # Use vegan package for equal-depth resampling, input reads count format Feature table result/otutab.txt.
    # Can specify input file, sampling amount, and random number, output rarefied table result/otutab_rare.txt and diversity alpha/vegan.txt.
    mkdir -p "${OUTPUT_DIR}/alpha"
    cp -f "${INPUT_OTUTAB_FILE}" "${TEMP_INPUT_OTUTAB_FILE}"
    Rscript "${db}/script/otutab_rare.R" --input "${TEMP_INPUT_OTUTAB_FILE}" \
      --depth "${DEPTH}" --seed 1 \
      --normalize "${OUTPUT_RAREFIED_OTUTAB_FILE}" \
      --output "${OUTPUT_ALPHA_DIVERSITY_FILE}"

    # Rarefied table stats (optional)
    if command -v usearch >/dev/null 2>&1 && usearch -version >/dev/null 2>&1; then
        usearch -otutab_stats "${OUTPUT_RAREFIED_OTUTAB_FILE}" -output "${OUTPUT_RAREFIED_OTUTAB_FILE}.stat"
        cat "${OUTPUT_RAREFIED_OTUTAB_FILE}.stat"
    else
        echo "Warning: usearch unavailable, skip otutab_stats for rarefied table." >&2
    fi
