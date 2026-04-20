#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.
# Module: diversity_alpha - Alpha Diversity Analysis
# Function: Calculates and visualizes alpha diversity indices.

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
DBDIR="/d/EasyMicrobiome"
THREADS=4
OUTPUT_DIR="result"
RAREFIED_OTU_TABLE_FILE="${OUTPUT_DIR}/otutab_rare.txt"
FEATURE_SEQS_FILE="${OUTPUT_DIR}/otus.fa" # New parameter
METADATA_FILE="${OUTPUT_DIR}/metadata.txt"
ALPHA_DIVERSITY_FILE="${OUTPUT_DIR}/alpha/alpha.txt" # Changed from alpha.txt7
ALPHA_RARE_FILE="${OUTPUT_DIR}/alpha/alpha_rare.txt"

# --- Help Documentation ---
usage() {
    cat <<EOF
Module: diversity_alpha - Alpha Diversity Analysis

Usage: bash easyamplicon.sh diversity_alpha [OPTIONS]

Options:
  -d, --dbdir DIR        Set the EasyMicrobiome database directory.
                         (Default: ${DBDIR})
  -t, --threads INT      Set the number of threads for parallel tasks.
                         (Default: ${THREADS})
  --output-dir DIR       Main results output directory.
                         (Default: ${OUTPUT_DIR})
  --rarefied-otu-table FILE Path to the rarefied OTU table file.
                         (Default: ${RAREFIED_OTU_TABLE_FILE})
  --feature-seqs FILE    Path to the feature sequences file (e.g., result/otus.fa).
                         (Default: ${FEATURE_SEQS_FILE})
  --metadata-file FILE   Path to the metadata file, containing sample information.
                         (Default: ${METADATA_FILE})
  --alpha-diversity-output FILE Path for the alpha diversity indices output.
                             (Default: ${ALPHA_DIVERSITY_FILE})
  --alpha-rare-output FILE Path for the rarefaction curve data output.
                             (Default: ${ALPHA_RARE_FILE})
  -h, --help             Display this help message and exit.
EOF
    exit 0
}

# --- Parse Command Line Arguments ---
TEMP=$(getopt -o d:t:h --long dbdir:,threads:,output-dir:,rarefied-otu-table:,feature-seqs:,metadata-file:,alpha-diversity-output:,alpha-rare-output:,help -n "$0" -- "$@")
if [ $? != 0 ]; then echo "Terminating..." >&2; exit 1; fi
eval set -- "$TEMP"

while true; do
    case "$1" in
        -d|--dbdir) DBDIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --rarefied-otu-table) RAREFIED_OTU_TABLE_FILE="$2"; shift 2 ;;
        --feature-seqs) FEATURE_SEQS_FILE="$2"; shift 2 ;;
        --metadata-file) METADATA_FILE="$2"; shift 2 ;;
        --alpha-diversity-output) ALPHA_DIVERSITY_FILE="$2"; shift 2 ;;
        --alpha-rare-output) ALPHA_RARE_FILE="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) echo "Internal error!"; exit 1 ;;
    esac
done

# --- Export variables for sub-scripts (if needed) ---
export db=${DBDIR}
export threads=${THREADS}

# Ensure output directory for alpha diversity exists
mkdir -p "$(dirname "${ALPHA_DIVERSITY_FILE}")"

# --- Module Core Logic ---


### 6.1. Calculate alpha diversity

    # Calculate 14 alpha diversity indices using USEARCH (Chao1 is erroneous, do not use).
    # details in http://www.drive5.com/usearch/manual/alpha_metrics.html
    usearch -alpha_div "${RAREFIED_OTU_TABLE_FILE}" \
      -output "${ALPHA_DIVERSITY_FILE}"

### 6.2. Calculate rarefaction richness

    # Rarefaction curve: number of OTUs from 1%-100% sequences, sampling without replacement each time.
    # Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
    usearch -alpha_div_rare "${RAREFIED_OTU_TABLE_FILE}" \
      -output "${ALPHA_RARE_FILE}" \
      -method without_replacement
    # Preview results
    head -n2 "${ALPHA_RARE_FILE}"
    # Handling non-numeric "-" values due to low sample sequencing depth, see FAQ 8.
    sed -i "s/-/\t0.0/g" "${ALPHA_RARE_FILE}"
