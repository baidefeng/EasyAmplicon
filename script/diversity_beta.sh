#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.
# Module: diversity_beta - Beta Diversity Analysis
# Function: Calculates and visualizes beta diversity distances.

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
METADATA_FILE="${OUTPUT_DIR}/metadata.txt" # Changed to use OUTPUT_DIR
OTUS_FA_FILE="${OUTPUT_DIR}/otus.fa"
OTUS_TREE_FILE="${OUTPUT_DIR}/otus.tree"
BETA_OUTPUT_DIR="${OUTPUT_DIR}/beta" # New parameter for the beta output directory

# --- Help Documentation ---
usage() {
    cat <<EOF
Module: diversity_beta - Beta Diversity Analysis

Usage: bash easyamplicon.sh diversity_beta [OPTIONS]

Options:
  -d, --dbdir DIR        Set the EasyMicrobiome database directory.
                         (Default: ${DBDIR})
  -t, --threads INT      Set the number of threads for parallel tasks.
                         (Default: ${THREADS})
  --output-dir DIR       Main results output directory.
                         (Default: ${OUTPUT_DIR})
  --rarefied-otu-table FILE Path to the rarefied OTU table file.
                         (Default: ${RAREFIED_OTU_TABLE_FILE})
  --metadata-file FILE   Path to the metadata file, containing sample information.
                         (Default: ${METADATA_FILE})
  --otus-fa FILE         Path to the feature sequences (OTUs/ASVs) FASTA file.
                         (Default: ${OTUS_FA_FILE})
  --otus-tree FILE       Path for the output phylogenetic tree file.
                         (Default: ${OTUS_TREE_FILE})
  --beta-output-dir DIR  Directory for beta diversity output files.
                         (Default: ${BETA_OUTPUT_DIR})
  -h, --help             Display this help message and exit.
EOF
    exit 0
}

# --- Parse Command Line Arguments ---
TEMP=$(getopt -o d:t:h --long dbdir:,threads:,output-dir:,rarefied-otu-table:,metadata-file:,otus-fa:,otus-tree:,beta-output-dir:,help -n "$0" -- "$@")
if [ $? != 0 ]; then echo "Terminating..." >&2; exit 1; fi
eval set -- "$TEMP"

while true; do
    case "$1" in
        -d|--dbdir) DBDIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --rarefied-otu-table) RAREFIED_OTU_TABLE_FILE="$2"; shift 2 ;;
        --metadata-file) METADATA_FILE="$2"; shift 2 ;;
        --otus-fa) OTUS_FA_FILE="$2"; shift 2 ;;
        --otus-tree) OTUS_TREE_FILE="$2"; shift 2 ;;
        --beta-output-dir) BETA_OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) echo "Internal error!"; exit 1 ;;
    esac
done

# --- Export variables for sub-scripts (if needed) ---
export db=${DBDIR}
export threads=${THREADS}

# --- Module Core Logic ---

## Beta Diversity Analysis

    echo "--- Beta Diversity Analysis ---"
    # Results have multiple files, directory is needed.
    mkdir -p "${BETA_OUTPUT_DIR}"
    # Make OTU tree, 4s
    usearch -cluster_agg "${OTUS_FA_FILE}" -treeout "${OTUS_TREE_FILE}"
    # Generate 5 distance matrices: bray_curtis, euclidean, jaccard, manhatten, unifrac
    usearch -beta_div "${RAREFIED_OTU_TABLE_FILE}" -tree "${OTUS_TREE_FILE}" \
      -filename_prefix "${BETA_OUTPUT_DIR}/"
