#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.
# Module: feature_table - Core Sequence Processing, OTU/ASV Generation, and Feature Table Construction
# Function: Integrates sequence merging, renaming, quality control, dereplication,
#           OTU/ASV generation, chimera removal, and feature table construction and filtering.

# --- Default Parameters ---
# EasyMicrobiome root (contains script/, usearch/, etc.). Override per machine:
#   export EASY_MICROBIOME=/path/to/EasyMicrobiome
#   or: bash easyamplicon.sh feature_table --dbdir /path/to/EasyMicrobiome ...
DBDIR="${EASY_MICROBIOME:-/d/EasyMicrobiome}"
THREADS=4

# Input/Output Paths
RAW_READS_DIR="seq"
METADATA_FILE="result/metadata.txt"
TEMP_DIR="temp"
OUTPUT_DIR="result"
OTUTAB_FILENAME="otutab.txt" # New parameter for otutab filename
OTUS_FILENAME="otus.fa"     # New parameter for otus filename
MERGED_READS_DIR="${TEMP_DIR}"
CONCAT_MERGED_READS_FILE="${TEMP_DIR}/all.fq"
FILTERED_READS_FILE="${TEMP_DIR}/filtered.fa"
DEREPLICATED_SEQS_FILE="${TEMP_DIR}/uniques.fa"
ASV_OTU_SEQS_FILE="${TEMP_DIR}/otus.fa"
CHIMERA_FILTERED_SEQS_FILE="${OUTPUT_DIR}/raw/otus.fa"
RAW_OTU_TABLE_FILE="${OUTPUT_DIR}/raw/otutab.txt"
RAW_TAXONOMY_FILE="${OUTPUT_DIR}/raw/otus.sintax"
FILTERED_OTU_TABLE_FILE="${OUTPUT_DIR}/${OTUTAB_FILENAME}"
FILTERED_OTU_SEQS_FILE="${OUTPUT_DIR}/${OTUS_FILENAME}"
FILTERED_TAXONOMY_FILE="${OUTPUT_DIR}/otus.sintax"
OTU_TABLE_STATS_FILE="${OUTPUT_DIR}/otutab.stat"
ALPHA_DIVERSITY_FILE="${OUTPUT_DIR}/alpha/vegan.txt"
RAREFIED_OTU_TABLE_FILE="${OUTPUT_DIR}/otutab_rare.txt"

# Processing Parameters
LEFT_STRIP=29
RIGHT_STRIP=18
MAXEE_RATE=0.01
MIN_UNIQUE_SIZE=10
CHIMERA_REF_DB="${DBDIR}/usearch/rdp_16s_v18.fa"
SINTAX_DB="${DBDIR}/usearch/rdp_16s_v18.fa"
SINTAX_CUTOFF=0.1
ID_CUTOFF=0.97
FILTER_NONBAC="TRUE" # Can be "TRUE" or "FALSE"
FILTER_NONBAC_R="${DBDIR}/script/otutab_filter_nonBac.R"

# --- Help Documentation ---
usage() {
    cat <<EOF
Module: feature_table - Core Sequence Processing, OTU/ASV Generation, and Feature Table Construction

Usage: bash easyamplicon.sh feature_table [OPTIONS]

Options:
  -d, --dbdir DIR        Set the EasyMicrobiome database directory.
                         (Default: ${DBDIR})
  -t, --threads INT      Set the number of threads for parallel tasks.
                         (Default: ${THREADS})

  --raw-reads-dir DIR    Directory containing raw FASTQ.GZ files.
                         (Default: ${RAW_READS_DIR})
  --metadata-file FILE   Path to the metadata file, containing sample information.
                         (Default: ${METADATA_FILE})
  --temp-dir DIR         Directory for temporary files.
                         (Default: ${TEMP_DIR})
  --output-dir DIR       Main results output directory.
                         (Default: ${OUTPUT_DIR})
  --otutab-filename FILE Name for the filtered feature table (OTU/ASV table).
                         (Default: ${OTUTAB_FILENAME})
  --otus-filename FILE   Name for the filtered feature sequences.
                         (Default: ${OTUS_FILENAME})

  --left-strip INT       Number of bases to strip from the left end of sequences (primer/barcode).
                         (Default: ${LEFT_STRIP})
  --right-strip INT      Number of bases to strip from the right end of sequences (primer).
                         (Default: ${RIGHT_STRIP})
  --maxee-rate FLOAT     Maximum expected error rate for quality filtering.
                         (Default: ${MAXEE_RATE})
  --min-unique-size INT  Minimum size for unique sequences during dereplication.
                         (Default: ${MIN_UNIQUE_SIZE})
  --chimera-ref-db FILE  Path to the reference database for chimera detection.
                         (Default: ${CHIMERA_REF_DB})
  --sintax-db FILE       SINTAX database for taxonomy annotation.
                         (Default: ${SINTAX_DB})
  --sintax-cutoff FLOAT  SINTAX confidence cutoff.
                         (Default: ${SINTAX_CUTOFF})
  --id-cutoff FLOAT      Similarity threshold for feature table generation.
                         (Default: ${ID_CUTOFF})
  --filter-nonbac BOOL   Filter non-bacteria/archaea, chloroplast, mitochondria (TRUE/FALSE).
                         (Default: ${FILTER_NONBAC})

  -h, --help             Display this help message and exit.
EOF
    exit 0
}

# --- Parse Command Line Arguments ---
TEMP=$(getopt -o d:t:h --long dbdir:,threads:,raw-reads-dir:,metadata-file:,temp-dir:,output-dir:,otutab-filename:,otus-filename:,left-strip:,right-strip:,maxee-rate:,min-unique-size:,chimera-ref-db:,sintax-db:,sintax-cutoff:,id-cutoff:,filter-nonbac:,help -n "$0" -- "$@")
if [ $? != 0 ]; then echo "Terminating..." >&2; exit 1; fi
eval set -- "$TEMP"

while true; do
    case "$1" in
        -d|--dbdir) DBDIR="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --raw-reads-dir) RAW_READS_DIR="$2"; shift 2 ;;
        --metadata-file) METADATA_FILE="$2"; shift 2 ;;
        --temp-dir) TEMP_DIR="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --otutab-filename) OTUTAB_FILENAME="$2"; shift 2 ;;
        --otus-filename) OTUS_FILENAME="$2"; shift 2 ;;
        --left-strip) LEFT_STRIP="$2"; shift 2 ;;
        --right-strip) RIGHT_STRIP="$2"; shift 2 ;;
        --maxee-rate) MAXEE_RATE="$2"; shift 2 ;;
        --min-unique-size) MIN_UNIQUE_SIZE="$2"; shift 2 ;;
        --chimera-ref-db) CHIMERA_REF_DB="$2"; shift 2 ;;
        --sintax-db) SINTAX_DB="$2"; shift 2 ;;
        --sintax-cutoff) SINTAX_CUTOFF="$2"; shift 2 ;;
        --id-cutoff) ID_CUTOFF="$2"; shift 2 ;;
        --filter-nonbac) FILTER_NONBAC="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) echo "Internal error!"; exit 1 ;;
    esac
done

# --- Export variables for sub-scripts (if needed) ---
export db=${DBDIR}
export threads=${THREADS}

# Resolve dependent paths after parsing args
FILTER_NONBAC_R="${DBDIR}/script/otutab_filter_nonBac.R"

# Basic runtime checks (keep light; downstream tools may live in conda env)
for cmd in vsearch Rscript csvtk awk cut sed; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "Error: '${cmd}' not found in PATH." >&2
        exit 1
    fi
done

# --- Helpers ---
require_file() {
    local f="$1"
    local what="${2:-file}"
    if [ ! -f "$f" ]; then
        echo "Error: missing ${what}: $f" >&2
        echo "  Download it yourself and place it at this path." >&2
        echo "  Download link (TODO): <PUT_DOWNLOAD_LINK_HERE>" >&2
        exit 1
    fi
}

require_nonempty() {
    local f="$1"
    local what="${2:-file}"
    if [ ! -s "$f" ]; then
        echo "Error: ${what} is missing or empty: $f" >&2
        exit 1
    fi
}

require_glob_any() {
    local pattern="$1"
    local what="${2:-files}"
    if ! ls ${pattern} >/dev/null 2>&1; then
        echo "Error: expected ${what} not found: ${pattern}" >&2
        exit 1
    fi
}

ensure_uncompressed() {
    # If target exists, do nothing.
    # If only target.gz exists, decompress to target (keep .gz).
    local target="$1"
    local what="${2:-file}"
    if [ -f "$target" ]; then
        return 0
    fi
    if [ -f "${target}.gz" ]; then
        if ! command -v gzip >/dev/null 2>&1; then
            echo "Error: gzip not found; cannot decompress ${target}.gz" >&2
            exit 1
        fi
        echo "Decompressing ${target}.gz -> ${target}"
        gzip -dc "${target}.gz" > "${target}"
        return 0
    fi
    echo "Error: missing ${what}: ${target} (or ${target}.gz)" >&2
    echo "  Download it yourself and place it at this path." >&2
    echo "  Download link (TODO): <PUT_DOWNLOAD_LINK_HERE>" >&2
    exit 1
}
if [ "${FILTER_NONBAC}" = "TRUE" ] && [ ! -f "${FILTER_NONBAC_R}" ]; then
    echo "Error: non-bacteria filter script not found: ${FILTER_NONBAC_R}" >&2
    echo "  Use --dbdir /path/to/EasyMicrobiome  or  export EASY_MICROBIOME=..." >&2
    echo "  Or  --filter-nonbac FALSE" >&2
    exit 1
fi

# --- Ensure key databases exist (supports .gz) ---
ensure_uncompressed "${CHIMERA_REF_DB}" "chimera reference database"
ensure_uncompressed "${SINTAX_DB}" "SINTAX taxonomy database"

# Update paths based on parameters
MERGED_READS_DIR="${TEMP_DIR}"
CONCAT_MERGED_READS_FILE="${TEMP_DIR}/all.fq"
FILTERED_READS_FILE="${TEMP_DIR}/filtered.fa"
DEREPLICATED_SEQS_FILE="${TEMP_DIR}/uniques.fa"
ASV_OTU_SEQS_FILE="${TEMP_DIR}/otus.fa"
CHIMERA_FILTERED_SEQS_FILE="${OUTPUT_DIR}/raw/otus.fa"
RAW_OTU_TABLE_FILE="${OUTPUT_DIR}/raw/otutab.txt"
RAW_TAXONOMY_FILE="${OUTPUT_DIR}/raw/otus.sintax"
FILTERED_OTU_TABLE_FILE="${OUTPUT_DIR}/${OTUTAB_FILENAME}"
FILTERED_OTU_SEQS_FILE="${OUTPUT_DIR}/${OTUS_FILENAME}"
FILTERED_TAXONOMY_FILE="${OUTPUT_DIR}/otus.sintax"
OTU_TABLE_STATS_FILE="${OUTPUT_DIR}/otutab.stat"
ALPHA_DIVERSITY_FILE="${OUTPUT_DIR}/alpha/vegan.txt"
RAREFIED_OTU_TABLE_FILE="${OUTPUT_DIR}/otutab_rare.txt"


# --- Module Core Logic ---

# Content from reads_processing.sh
# Module 2: Sequence Merging, Renaming, and Quality Control
# Function: Merges paired-end sequences, renames by sample, concatenates sequences,
#           and performs primer trimming and quality filtering.

## Reads merge and rename

### Merge paired-end reads and rename by sample


    # Process sequentially using a for loop (if rush is not available)
    echo "Processing sequentially using a for loop (if rush is not available)"
    require_file "${METADATA_FILE}" "metadata file"
    time for i in `tail -n+2 "${METADATA_FILE}"|cut -f1`;do
      require_file "${RAW_READS_DIR}/${i}_1.fq.gz" "raw reads (forward) for sample ${i}"
      require_file "${RAW_READS_DIR}/${i}_2.fq.gz" "raw reads (reverse) for sample ${i}"
      vsearch --fastq_mergepairs "${RAW_READS_DIR}/${i}_1.fq.gz" --reverse "${RAW_READS_DIR}/${i}_2.fq.gz" \
      --fastqout "${MERGED_READS_DIR}/${i}.merged.fq" --relabel "${i}."
    done


### (Optional) Single-end reads rename

    # # Example of renaming a single sequence
    # i=WT1
    # gunzip -c "${RAW_READS_DIR}/${i}_1.fq.gz" > "${RAW_READS_DIR}/${i}.fq"
    # usearch -fastx_relabel "${RAW_READS_DIR}/${i}.fq" -fastqout "${MERGED_READS_DIR}/${i}.merged.fq" -prefix "${i}."
    # 
    # # Batch rename, requires single-end fastq files, and decompressed (usearch does not support compressed format)
    # gunzip "${RAW_READS_DIR}"/*.gz
    # time for i in `tail -n+2 "${METADATA_FILE}"|cut -f1`;do
    #   usearch -fastx_relabel "${RAW_READS_DIR}/${i}.fq" -fastqout "${MERGED_READS_DIR}/${i}.merged.fq" -prefix "${i}."
    # done &
    # # For large datasets, refer to "FAQ 2" for the vsearch method.

### Integrate renamed reads

    echo "--- Integrate renamed reads ---"
    if ! ls "${MERGED_READS_DIR}"/*.merged.fq >/dev/null 2>&1; then
        echo "Error: no merged reads found at ${MERGED_READS_DIR}/*.merged.fq" >&2
        exit 1
    fi
    # Merge all samples into a single file
    cat "${MERGED_READS_DIR}"/*.merged.fq > "${CONCAT_MERGED_READS_FILE}"
    require_nonempty "${CONCAT_MERGED_READS_FILE}" "integrated merged FASTQ"
    # Check the file size (e.g., 223M). Results may vary slightly with different software versions.
    ls -lsh "${CONCAT_MERGED_READS_FILE}"
    # Check the sequence names. The part before the "." should be the sample name. Sample names must not contain dots (".").
    # A significant feature of sample names with dots is that the generated feature table will be very large, with many columns, leading to memory shortages in subsequent analysis.
    # After obtaining the feature table, you should check it for any issues. If you encounter memory problems, you should go back and investigate.
    head -n 6 "${CONCAT_MERGED_READS_FILE}"|cut -c1-60


## Cut primers and quality filter

    echo "--- Cut primers and quality filter ---"
    require_file "${CONCAT_MERGED_READS_FILE}" "merged reads file"
    # 10bp barcode + 19bp upstream primer V5 total 29 from left, 18bp downstream primer V7 from right.
    # Be sure to understand the experimental design and primer lengths. If primers have already been removed, set to 0. Takes 14s for 270k sequences.
    time vsearch --fastx_filter "${CONCAT_MERGED_READS_FILE}" \
      --fastq_stripleft "${LEFT_STRIP}" --fastq_stripright "${RIGHT_STRIP}" \
      --fastq_maxee_rate "${MAXEE_RATE}" \
      --fastaout "${FILTERED_READS_FILE}"
    require_nonempty "${FILTERED_READS_FILE}" "filtered reads FASTA"
    # View the file to understand the FASTA format
    head "${FILTERED_READS_FILE}"

# Content from otu_asv_generation.sh
# Module: otu_asv_generation - Dereplicate and Select OTU/ASV
# Function: Dereplicates, clusters or denoises quality-controlled sequences, and removes chimeras to generate representative sequences.

## Dereplicate and cluster/denoise

### Dereplicate sequences

    echo "--- Dereplicate sequences ---"
    require_file "${FILTERED_READS_FILE}" "filtered reads fasta"
    # Add minuniquesize (minimum 8 or 1/1M) to remove low-abundance noise and increase calculation speed.
    # -sizeout outputs abundance, --relabel must add a sequence prefix for better standardization. Takes 1s.
    vsearch --derep_fulllength "${FILTERED_READS_FILE}" \
      --minuniquesize "${MIN_UNIQUE_SIZE}" --sizeout --relabel Uni_ \
      --output "${DEREPLICATED_SEQS_FILE}"
    require_nonempty "${DEREPLICATED_SEQS_FILE}" "dereplicated sequences FASTA"
    # High-abundance non-redundant sequences are very small (500K~5M is more suitable), with size and frequency after the name.
    ls -lsh "${DEREPLICATED_SEQS_FILE}"
    # Uni_1;size=6423 - The name of the sequence after dereplication is Uni_1; this sequence appears 6423 times in all sample sequencing data.
    # This is the most abundant sequence.
    head -n 2 "${DEREPLICATED_SEQS_FILE}"

### Cluster OTUs / Denoise ASVs

    echo "--- Cluster OTUs / Denoise ASVs ---"
    require_file "${DEREPLICATED_SEQS_FILE}" "dereplicated sequences fasta"
    # There are two methods: unoise3 denoising to obtain single-base precision ASVs is recommended. The alternative is traditional 97% clustering for OTUs (genus-level precision).
    # Both feature selection methods in usearch include de novo chimera removal.
    # -minsize for secondary filtering, to control the number of OTUs/ASVs to 1-5 thousand, which is convenient for downstream statistical analysis.

    # Method 1. 97% OTU clustering, suitable for big data / when ASV patterns are not obvious / required by reviewers.
    # Results took 1s, produced 508 OTUs, and removed 126 chimeras.
    # usearch -cluster_otus "${DEREPLICATED_SEQS_FILE}" -minsize "${MIN_UNIQUE_SIZE}" \
    #  -otus "${ASV_OTU_SEQS_FILE}" \
    #  -relabel OTU_

    # Method 2. ASV Denoise: predict biological sequences and filter chimeras
    # Prefer usearch UNOISE3 when available; otherwise fallback to vsearch cluster_unoise.
    HAS_USEARCH="TRUE"
    if ! command -v usearch >/dev/null 2>&1; then
        HAS_USEARCH="FALSE"
    elif ! usearch -version >/dev/null 2>&1; then
        HAS_USEARCH="FALSE"
    fi

    if [ "${HAS_USEARCH}" = "TRUE" ]; then
        usearch -unoise3 "${DEREPLICATED_SEQS_FILE}" -minsize "${MIN_UNIQUE_SIZE}" \
          -zotus "${TEMP_DIR}/zotus.fa"
    else
        echo "Warning: usearch unavailable, fallback to vsearch --cluster_unoise." >&2
        vsearch --cluster_unoise "${DEREPLICATED_SEQS_FILE}" \
          --minsize "${MIN_UNIQUE_SIZE}" \
          --centroids "${TEMP_DIR}/zotus.fa"
    fi

    # Normalize IDs:
    # - Convert Uni_/Centroid_/Zotu_ to ASV_
    # - Strip ";size=..." to make IDs match feature table rownames and downstream R filters
    sed -E '
      s/^>Zotu(_)?/>ASV_/;
      s/^>Centroid_/>ASV_/;
      s/^>Uni_/>ASV_/;
      s/^>([^;[:space:]]+);size=[0-9]+.*/>\1/;
    ' "${TEMP_DIR}/zotus.fa" > "${ASV_OTU_SEQS_FILE}"
    require_nonempty "${ASV_OTU_SEQS_FILE}" "ASV representative sequences FASTA"
    head -n 2 "${ASV_OTU_SEQS_FILE}"

    # Method 3. When the data is too large to use usearch, see "FAQ 3" for the alternative vsearch method.

### Reference-based chimera detection

    echo "--- Reference-based chimera detection ---"
    # Not recommended, as it can easily cause false negatives because the reference database lacks abundance information.
    # In de novo chimera detection, the abundance of the parent sequences is required to be more than 16 times that of the chimera to prevent false negatives.
    # Since known sequences will not be removed, the larger the database selected, the more reasonable it is, and the lower the false negative rate.
    mkdir -p "${OUTPUT_DIR}/raw"

    # Method 1. Chimera removal with vsearch + rdp (fast but prone to false negatives)
    vsearch --uchime_ref "${ASV_OTU_SEQS_FILE}" \
      -db "${CHIMERA_REF_DB}" \
      --nonchimeras "${CHIMERA_FILTERED_SEQS_FILE}"
    require_nonempty "${CHIMERA_FILTERED_SEQS_FILE}" "chimera-filtered feature FASTA"
    # RDP: 7s, 143 (9.3%) chimeras; SILVA：9m, 151 (8.4%) chimeras
    # The results of vsearch on Windows have added a Windows newline character (^M), which needs to be deleted. Do not execute this command on a Mac.
    sed -i 's/\r//g' "${CHIMERA_FILTERED_SEQS_FILE}"

    # Method 2. Do not remove chimeras
    # cp -f "${ASV_OTU_SEQS_FILE}" "${CHIMERA_FILTERED_SEQS_FILE}"

# Content from feature_table.sh
# Module: feature_table - Feature Table Construction and Filtering
# Function: Generates feature table, performs taxonomic annotation, removes contaminant sequences, and performs sampling normalization.

## Feature table create and filter

# OTUs and ASVs are collectively referred to as Features. Their differences are:
# OTUs are usually representative sequences selected with the highest abundance or from the center after 97% clustering.
# ASVs are representative sequences obtained by denoising based on sequences (excluding or correcting erroneous sequences and selecting credible sequences with higher abundance).

### Generate a Feature table

    echo "--- Generate a Feature table ---"
    require_file "${FILTERED_READS_FILE}" "filtered reads fasta"
    require_file "${CHIMERA_FILTERED_SEQS_FILE}" "chimera-filtered feature fasta"
    # id(1): 100% similarity alignment, 49.45% sequences, 1m50s
    # id(0.97): 97% similarity alignment, 83.66% sequences, 1m10s (higher data usage, faster)
    time vsearch --usearch_global "${FILTERED_READS_FILE}" \
      --db "${CHIMERA_FILTERED_SEQS_FILE}" \
      --id "${ID_CUTOFF}" --threads "${THREADS}" \
    	--otutabout "${RAW_OTU_TABLE_FILE}"
    require_nonempty "${RAW_OTU_TABLE_FILE}" "raw feature table"
    # 212862 of 268019 (79.42%) can be aligned
    # For vsearch results, Windows users should delete the newline character (^M) to correct it to the standard Linux format.
    sed -i 's/\r//' "${RAW_OTU_TABLE_FILE}"
    head -n6 "${RAW_OTU_TABLE_FILE}" | cut -f 1-6 |cat -A
    # Use csvtk to count the rows and columns of the table.
    # Be sure to check the number of columns here. Is it equal to your number of samples? If not, there is generally a problem with the sample naming. See the explanation above for details.
    csvtk -t stat "${RAW_OTU_TABLE_FILE}"

    # Normalize OTU IDs by stripping ";size=..." so that:
    # - Feature table rownames (ASV_1) match
    # - SINTAX first column and downstream R filters match
    CHIMERA_FILTERED_SEQS_NOSIZE_FILE="${OUTPUT_DIR}/raw/otus_nosize.fa"
    sed -E 's/^>([^;[:space:]]+);size=[0-9]+.*/>\1/' "${CHIMERA_FILTERED_SEQS_FILE}" > "${CHIMERA_FILTERED_SEQS_NOSIZE_FILE}"
    require_nonempty "${CHIMERA_FILTERED_SEQS_NOSIZE_FILE}" "raw feature sequences (nosize)"

    echo "--- Taxonomic Annotation and Removal of Plastids and Non-Bacteria ---"
    require_file "${RAW_OTU_TABLE_FILE}" "raw feature table"
    require_file "${CHIMERA_FILTERED_SEQS_NOSIZE_FILE}" "raw feature sequences (nosize)"
    vsearch --sintax "${CHIMERA_FILTERED_SEQS_NOSIZE_FILE}" \
      --db "${SINTAX_DB}" \
      --sintax_cutoff "${SINTAX_CUTOFF}" \
      --tabbedout "${RAW_TAXONOMY_FILE}"
    sed -i 's/\r//g' "${RAW_TAXONOMY_FILE}"
    require_nonempty "${RAW_TAXONOMY_FILE}" "raw taxonomy table"

    if [ "${FILTER_NONBAC}" = "TRUE" ]; then
        Rscript "${FILTER_NONBAC_R}" \
          --input "${RAW_OTU_TABLE_FILE}" \
          --taxonomy "${RAW_TAXONOMY_FILE}" \
          --output "${FILTERED_OTU_TABLE_FILE}" \
          --stat "${OUTPUT_DIR}/raw/otutab_nonBac.stat" \
          --discard "${OUTPUT_DIR}/raw/otus.sintax.discard"
        require_nonempty "${FILTERED_OTU_TABLE_FILE}" "filtered feature table"

        # Filter sequences to match filtered table
        cut -f 1 "${FILTERED_OTU_TABLE_FILE}" | tail -n+2 > "${OUTPUT_DIR}/otutab.id"
        require_nonempty "${OUTPUT_DIR}/otutab.id" "filtered feature ID list"
        vsearch --fastx_getseqs "${CHIMERA_FILTERED_SEQS_NOSIZE_FILE}" \
          --labels "${OUTPUT_DIR}/otutab.id" \
          --fastaout "${FILTERED_OTU_SEQS_FILE}"
        require_nonempty "${FILTERED_OTU_SEQS_FILE}" "filtered feature sequences FASTA"

        # Filter taxonomy to match filtered table
        awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}' \
          "${RAW_TAXONOMY_FILE}" "${OUTPUT_DIR}/otutab.id" \
          > "${FILTERED_TAXONOMY_FILE}"
        require_nonempty "${FILTERED_TAXONOMY_FILE}" "filtered taxonomy table"
    else
        # No filtering: just publish normalized raw outputs
        cp -f "${RAW_OTU_TABLE_FILE}" "${FILTERED_OTU_TABLE_FILE}"
        cp -f "${CHIMERA_FILTERED_SEQS_NOSIZE_FILE}" "${FILTERED_OTU_SEQS_FILE}"
        cp -f "${RAW_TAXONOMY_FILE}" "${FILTERED_TAXONOMY_FILE}"
        require_nonempty "${FILTERED_OTU_TABLE_FILE}" "feature table"
        require_nonempty "${FILTERED_OTU_SEQS_FILE}" "feature sequences FASTA"
        require_nonempty "${FILTERED_TAXONOMY_FILE}" "taxonomy table"
    fi

    # Optional summary stats (usearch provides otutab_stats; if missing, skip)
    if command -v usearch >/dev/null 2>&1 && usearch -version >/dev/null 2>&1; then
        usearch -otutab_stats "${FILTERED_OTU_TABLE_FILE}" -output "${OTU_TABLE_STATS_FILE}"
    else
        echo "Warning: usearch unavailable, skip otutab_stats (${OTU_TABLE_STATS_FILE})." >&2
    fi
    
    echo "--- [feature_table] Done ---"