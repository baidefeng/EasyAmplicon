#!/bin/bash

# Module: compare - Differential Comparison Analysis
# Function: Performs differential abundance analysis using R, STAMP, and LEfSe,
#           and generates corresponding input files and visualization results.

# --- Default Parameters ---
DBDIR="/d/EasyMicrobiome"
GROUP="Group"
COMPARE="KO-WT"
OUTPUT_DIR="result" # General output directory

FEATURE_TABLE_FILE="${OUTPUT_DIR}/otutab.txt"
METADATA_FILE="${OUTPUT_DIR}/metadata.txt"
TAXONOMY_FILE="${OUTPUT_DIR}/taxonomy.txt"
TAX_SUMMARY_DIR="${OUTPUT_DIR}/tax" # Base directory for sum_*.txt files

COMPARE_OUTPUT_DIR="${OUTPUT_DIR}/compare"
STAMP_OUTPUT_DIR="${OUTPUT_DIR}/stamp"
LEFSE_OUTPUT_DIR="${OUTPUT_DIR}/lefse"

# --- Help Documentation ---
usage() {
    cat <<EOF
Module: compare - Differential Comparison Analysis

Usage: bash easyamplicon.sh compare [OPTIONS]

Options:
  -d, --dbdir DIR        Set the EasyMicrobiome database directory.
                         (Default: /d/EasyMicrobiome)
  --group STRING         Set the group column name in the metadata.
                         (Default: "Group")
  --compare STRING       Set the comparison groups for differential analysis (e.g., "KO-WT").
                         (Default: "KO-WT")
  --feature-table FILE   Path to the filtered feature table file.
                         (Default: ${FEATURE_TABLE_FILE})
  --metadata-file FILE   Path to the metadata file.
                         (Default: ${METADATA_FILE})
  --taxonomy-file FILE   Path to the taxonomic annotation file.
                         (Default: ${TAXONOMY_FILE})
  --tax-summary-dir DIR  Directory containing summarized taxonomic abundance tables (sum_*.txt).
                         (Default: ${TAX_SUMMARY_DIR})
  --compare-output-dir DIR Directory for differential comparison analysis results.
                         (Default: ${COMPARE_OUTPUT_DIR})
  --stamp-output-dir DIR Directory for STAMP input and results.
                         (Default: ${STAMP_OUTPUT_DIR})
  --lefse-output-dir DIR Directory for LEfSe input files.
                         (Default: ${LEFSE_OUTPUT_DIR})
  -h, --help             Display this help message and exit.
EOF
    exit 0
}

# --- Parse Command Line Arguments ---
TEMP=$(getopt -o d:h --long dbdir:,group:,compare:,feature-table:,metadata-file:,taxonomy-file:,tax-summary-dir:,compare-output-dir:,stamp-output-dir:,lefse-output-dir:,help -n "$0" -- "$@")
if [ $? != 0 ]; then echo "Terminating..." >&2; exit 1; fi
eval set -- "$TEMP"

while true; do
    case "$1" in
        -d|--dbdir) DBDIR="$2"; shift 2 ;;
        --group) GROUP="$2"; shift 2 ;;
        --compare) COMPARE="$2"; shift 2 ;;
        --feature-table) FEATURE_TABLE_FILE="$2"; shift 2 ;;
        --metadata-file) METADATA_FILE="$2"; shift 2 ;;
        --taxonomy-file) TAXONOMY_FILE="$2"; shift 2 ;;
        --tax-summary-dir) TAX_SUMMARY_DIR="$2"; shift 2 ;;
        --compare-output-dir) COMPARE_OUTPUT_DIR="$2"; shift 2 ;;
        --stamp-output-dir) STAMP_OUTPUT_DIR="$2"; shift 2 ;;
        --lefse-output-dir) LEFSE_OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) echo "Internal error!"; exit 1 ;;
    esac
done

# --- Export variables for sub-scripts (if needed) ---
export db=${DBDIR}
export group=${GROUP}
export compare=${COMPARE}

# --- Module Core Logic ---

# 24. Differential Comparison

## 1. R-based Differential Analysis

### 1.1 Differential Comparison

    # Error in file(file, ifelse(append, "a", "w")), if output directory does not exist, create it.
    mkdir -p "${COMPARE_OUTPUT_DIR}"
    # Input feature table, metadata; specify group column name, comparison groups, and abundance.
    # Select method wilcox/t.test/edgeR, pvalue, fdr, and output directory.
Rscript ${db}/script/compare.R \
      --input "${FEATURE_TABLE_FILE}" --design "${METADATA_FILE}" \
      --group "${group}" --compare "${compare}" --threshold 0.1 \
      --method edgeR --pvalue 0.05 --fdr 0.2 \
      --output "${COMPARE_OUTPUT_DIR}/"

### 1.2 Volcano Plot

    # Input results from compare.R, output volcano plot with data labels, image size can be specified.
Rscript ${db}/script/compare_volcano.R \
      --input "${COMPARE_OUTPUT_DIR}/${compare}.txt" \
      --output "${COMPARE_OUTPUT_DIR}/${compare}.volcano.pdf" \
      --width 89 --height 59

### 1.3 Heatmap

    # Input results from compare.R, filter columns, specify metadata and groups, taxonomy annotation, image size in inches, and font size.
bash ${db}/script/compare_heatmap.sh -i "${COMPARE_OUTPUT_DIR}/${compare}.txt" -l 7 \
       -d "${METADATA_FILE}" -A "${group}" \
       -t "${TAXONOMY_FILE}" \
       -w 8 -h 5 -s 7 \
       -o "${COMPARE_OUTPUT_DIR}/${compare}"

### 1.4 Manhattan Plot

    # i: differential comparison results, t: taxonomy annotation, p: legend, w: width, v: height, s: font size, l: max legend value
    # If legend is not displayed, increase height v to 119+; later use AI to combine into KO-WT.heatmap.emf.
bash ${db}/script/compare_manhattan.sh -i "${COMPARE_OUTPUT_DIR}/${compare}.txt" \
       -t "${TAXONOMY_FILE}" \
       -p "${TAX_SUMMARY_DIR}/sum_p.txt" \
       -w 183 -v 59 -s 7 -l 10 \
       -o "${COMPARE_OUTPUT_DIR}/${compare}.manhattan.p.pdf"
    # The above plot only shows 6 phyla, switch to Class c and -L Class to show details.
bash ${db}/script/compare_manhattan.sh -i "${COMPARE_OUTPUT_DIR}/${compare}.txt" \
       -t "${TAXONOMY_FILE}" \
       -p "${TAX_SUMMARY_DIR}/sum_c.txt" \
       -w 183 -v 59 -s 7 -l 10 -L Class \
       -o "${COMPARE_OUTPUT_DIR}/${compare}.manhattan.c.pdf"
    # Display full legend, then use AI to combine.
bash ${db}/script/compare_manhattan.sh -i "${COMPARE_OUTPUT_DIR}/${compare}.txt" \
       -t "${TAXONOMY_FILE}" \
       -p "${TAX_SUMMARY_DIR}/sum_c.txt" \
       -w 183 -v 149 -s 7 -l 10 -L Class \
       -o "${COMPARE_OUTPUT_DIR}/${compare}.manhattan.c.legend.pdf"

### 1.5 Plotting Individual Features

    # Filter and display differential ASVs, sort by KO group abundance in descending order, take top 10 IDs.
awk '$4<0.05' "${COMPARE_OUTPUT_DIR}/${compare}.txt" | sort -k7,7nr | cut -f1 | head
    # Differential OTU details display
Rscript ${db}/script/alpha_boxplot.R --alpha_index ASV_2 \
      --input "${FEATURE_TABLE_FILE}" --design "${METADATA_FILE}" \
      --transpose TRUE --scale TRUE \
      --width 89 --height 59 \
      --group "${group}" --output "${COMPARE_OUTPUT_DIR}/feature_"
    # If ID does not exist, it will report an error: Error in data.frame(..., check.names = FALSE) : arguments imply differing number of rows: 0, 18  Calls: alpha_boxplot -> cbind -> cbind -> data.frame

    # Sort by a specific column: by mean abundance of genus All in descending order.
csvtk -t sort -k All:nr "${TAX_SUMMARY_DIR}/sum_g.txt" | head
    # Differential genus details display
Rscript ${db}/script/alpha_boxplot.R --alpha_index Lysobacter \
      --input "${TAX_SUMMARY_DIR}/sum_g.txt" --design "${METADATA_FILE}" \
      --transpose TRUE \
      --width 89 --height 59 \
      --group "${group}" --output "${COMPARE_OUTPUT_DIR}/feature_"

### 1.5 Ternary Plot

  # Refer to example in: result\compare\ternary\ternary.Rmd document
  # Alternative tutorial [246. Application and plotting practice of ternary plot](https://mp.weixin.qq.com/s/3w3ncpwjQaMRtmIOtr2Jvw)


## 2. STAMP Input File Preparation

### 2.1 Generate Input File

    Rscript ${db}/script/format2stamp.R -h
    mkdir -p "${STAMP_OUTPUT_DIR}"
Rscript ${db}/script/format2stamp.R --input "${FEATURE_TABLE_FILE}" \
      --taxonomy "${TAXONOMY_FILE}" --threshold 0.01 \
      --output "${STAMP_OUTPUT_DIR}/tax"
    # Optional Rmd document in result/format2stamp.Rmd

### 2.2 Plot Extended Bar Chart and Table

    # Replace ASV (result/otutab.txt) with genus (result/tax/sum_g.txt)
    # Code optimization: add legend, experiment vs control, bar order from top to bottom, currently WT top KO bottom.
Rscript ${db}/script/compare_stamp.R \
      --input "${STAMP_OUTPUT_DIR}/tax_5Family.txt" --metadata "${METADATA_FILE}" \
      --group "${group}" --compare "${compare}" --threshold 0.1 \
      --method "t.test" --pvalue 0.05 --fdr "none" \
      --width 280 --height 159 \
      --output "${STAMP_OUTPUT_DIR}/${compare}"
    # Optional Rmd document in result/CompareStamp.Rmd

## 3. LEfSe Input File Preparation

    ### 3.1. Generate File via Command Line
    # Optional command line to generate input file
    Rscript ${db}/script/format2lefse.R -h
    mkdir -p "${LEFSE_OUTPUT_DIR}"
    # threshold controls abundance filtering to control the number of branches in the plot.
Rscript ${db}/script/format2lefse.R --input "${FEATURE_TABLE_FILE}" \
      --taxonomy "${TAXONOMY_FILE}" --design "${METADATA_FILE}" \
      --group "${group}" --threshold 0.4 \
      --output "${LEFSE_OUTPUT_DIR}/LEfSe"

    ### 3.2 Generate Input File via Rmd (Optional)
    #1. otutab.txt, metadata.txt, taxonomy.txt three files exist in the result directory;
    #2. Open format2lefse.Rmd in EasyAmplicon with Rstudio, save to result directory and Knit to generate input file and reproducible computing webpage;

    ### 3.3 LEfSe Analysis
    #Method 1. Open LEfSe.txt and submit online to https://www.bic.ac.cn/BIC/#/analysis?tool_type=tool&page=b%27MzY%3D%27
    #Method 2. LEfSe local analysis (Linux system only, optional), refer to code in appendix
    #Method 3. LEfSe official website online use