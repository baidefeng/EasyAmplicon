#!/bin/bash

# EasyAmplicon v2.0 - Main Dispatcher Script

# --- Help Documentation ---
usage() {
    cat <<EOF
EasyAmplicon v2.0 - Modular Analysis Pipeline Dispatcher

Usage: $0 <module_name> [module_parameters]

Runs the specified analysis module and passes the corresponding parameters.

Available Modules:
  feature_table:      Core sequence processing, OTU/ASV generation, and feature table construction
  feature_pacbio:     PacBio feature table construction
  feature_nanopore:   Nanopore feature table construction
  feature_filter:     Feature table filtering and normalization
  taxonomy:           Taxonomic composition analysis
  diversity_alpha:    Alpha diversity analysis
  diversity_beta:     Beta diversity analysis
  compare:            Differential analysis (R, STAMP, LEfSe)
  function:           Functional prediction (PICRUSt, FAPROTAX, BugBase)


To view options for a specific module, use:
  bash $0 <module_name> --help
EOF
    exit 0
}

# --- Main Logic ---

# Display usage if no arguments or help is requested
if [ $# -eq 0 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    usage
fi

MODULE=$1
shift # Remove module name, keep its parameters

MODULE_SCRIPT="script/${MODULE}.sh"

# Check if the module script exists
if [[ -f "${MODULE_SCRIPT}" ]]; then
    # Strip CRLF without editing the file in place (avoids sed -i needing write
    # permission / temp files inside pipeline_modules on some mounts).
    bash <(sed 's/\r$//' "${MODULE_SCRIPT}") "$@"
else
    echo "Error: Unknown module '${MODULE}'"
    exit 1
fi
