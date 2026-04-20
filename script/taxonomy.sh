#!/bin/bash

# Module: taxonomy - Taxonomic Annotation and Classification Summary
# Function: Formats taxonomic annotations and summarizes them at different classification levels.

# --- Default Parameters ---
DBDIR="/d/EasyMicrobiome"
SINTAX_INPUT_FILE="result/otus.sintax"
TAX_OUTPUT_DIR="result/tax"
RAREFIED_OTU_TABLE_FILE="result/otutab_rare.txt"

# --- Help Documentation ---
usage() {
    cat <<EOF
Module: taxonomy - Taxonomic Annotation and Classification Summary

Usage: bash easyamplicon.sh taxonomy [OPTIONS]

Options:
  -d, --dbdir DIR        Set the EasyMicrobiome database directory.
                         (Default: /d/EasyMicrobiome)
  --sintax-input FILE    Path to the input sintax file (e.g., result/otus.sintax).
                         (Default: ${SINTAX_INPUT_FILE})
  --tax-output-dir DIR   Directory for taxonomic summary output (e.g., result/tax).
                         (Default: ${TAX_OUTPUT_DIR})
  --rarefied-otu-table FILE Path to the rarefied OTU table file (e.g., result/otutab_rare.txt).
                             (Default: ${RAREFIED_OTU_TABLE_FILE})
  -h, --help             Display this help message and exit.
EOF
    exit 0
}

# --- Parse Command Line Arguments ---
TEMP=$(getopt -o d:h --long dbdir:,sintax-input:,tax-output-dir:,rarefied-otu-table:,help -n "$0" -- "$@")
if [ $? != 0 ]; then echo "Terminating..." >&2; exit 1; fi
eval set -- "$TEMP"

while true; do
    case "$1" in
        -d|--dbdir) DBDIR="$2"; shift 2 ;;
        --sintax-input) SINTAX_INPUT_FILE="$2"; shift 2 ;;
        --tax-output-dir) TAX_OUTPUT_DIR="$2"; shift 2 ;;
        --rarefied-otu-table) RAREFIED_OTU_TABLE_FILE="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) echo "Internal error!"; exit 1 ;;
    esac
done

# --- Export variables for sub-scripts (if needed) ---
export db=${DBDIR}

# --- Module Core Logic ---

## 8. Taxonomic Annotation and Classification Summary

    # OTU corresponding to taxonomic annotation in 2-column format: remove confidence values from sintax,
    # only keep taxonomic annotation, replace ":" with "_", remove quotes.
    cut -f 1,4 "${SINTAX_INPUT_FILE}" \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    # OTU corresponding to taxonomic annotation in 8-column format: note that annotation may not be perfectly aligned.
    # Generate taxonomic table where blank OTU/ASV entries are filled with "Unassigned".
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt

    # Summarize Phylum, Class, Order, Family, Genus using rank parameters p c o f g,
    # which are abbreviations for phylum, class, order, family, genus.
    mkdir -p "${TAX_OUTPUT_DIR}"
    for i in p c o f g;do
      usearch -sintax_summary "${SINTAX_INPUT_FILE}" \
      -otutabin "${RAREFIED_OTU_TABLE_FILE}" -rank ${i} \
      -output "${TAX_OUTPUT_DIR}/sum_${i}.txt"
    done
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' "${TAX_OUTPUT_DIR}"/sum_*.txt
    # List all files
    wc -l "${TAX_OUTPUT_DIR}"/sum_*.txt
    head -n3 "${TAX_OUTPUT_DIR}"/sum_g.txt
