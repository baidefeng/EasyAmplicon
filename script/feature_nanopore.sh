#!/bin/bash
set -e
# Module: feature_nanopore - Nanopore feature table construction (stop before normalization)

# --- Default Parameters ---
DBDIR="${EASY_MICROBIOME:-/d/EasyMicrobiome}"
THREADS=4

RAW_READS_DIR="seq"
METADATA_FILE="result/metadata.txt"
TEMP_DIR="temp"
OUTPUT_DIR="result"
OTUTAB_FILENAME="otutab.txt"
OTUS_FILENAME="otus.fa"

PRIMER_FORWARD="AGAGTTTGATCCTGGCTCAG"
PRIMER_REVERSE="AAGTCSTAACAAGGTADCCSTA"
CUTADAPT_ERROR_RATE=0.1
MIN_LEN=1200
MAX_LEN=1800
FASTQ_MAXEE=30
FASTQ_QMAX=93
MIN_UNIQUE_SIZE=1
ID_CUTOFF=0.97
SINTAX_DB=""
SINTAX_CUTOFF=0.1
FILTER_NONBAC="TRUE"

# --- Help ---
usage() {
    cat <<EOF
Module: feature_nanopore - Nanopore feature table construction

Usage: bash easyamplicon.sh feature_nanopore [OPTIONS]

Options:
  -d, --dbdir DIR            EasyMicrobiome root. Default from \$EASY_MICROBIOME or /d/EasyMicrobiome.
                             (Current: ${DBDIR})
  -t, --threads INT          Thread count. (Default: ${THREADS})
  --raw-reads-dir DIR        Raw reads directory. (Default: ${RAW_READS_DIR})
  --metadata-file FILE       Metadata file path. (Default: ${METADATA_FILE})
  --temp-dir DIR             Temp directory. (Default: ${TEMP_DIR})
  --output-dir DIR           Output directory. (Default: ${OUTPUT_DIR})
  --otutab-filename FILE     Final feature table filename. (Default: ${OTUTAB_FILENAME})
  --otus-filename FILE       Final feature fasta filename. (Default: ${OTUS_FILENAME})
  --primer-forward STR       Forward primer. (Default: ${PRIMER_FORWARD})
  --primer-reverse STR       Reverse primer. (Default: ${PRIMER_REVERSE})
  --cutadapt-error-rate F    Cutadapt error rate. (Default: ${CUTADAPT_ERROR_RATE})
  --min-len INT              Minimum read length. (Default: ${MIN_LEN})
  --max-len INT              Maximum read length. (Default: ${MAX_LEN})
  --fastq-maxee FLOAT        vsearch max expected errors. (Default: ${FASTQ_MAXEE})
  --fastq-qmax INT           FASTQ quality max (ONT/PacBio: 93). (Default: ${FASTQ_QMAX})
  --min-unique-size INT      Minimum unique size for derep/cluster. (Default: ${MIN_UNIQUE_SIZE})
  --id-cutoff FLOAT          Mapping identity cutoff. (Default: ${ID_CUTOFF})
  --sintax-db FILE           SINTAX DB. (Default: ${SINTAX_DB})
  --sintax-cutoff FLOAT      SINTAX cutoff. (Default: ${SINTAX_CUTOFF})
  --filter-nonbac BOOL       Run non-bacteria filtering (TRUE/FALSE). (Default: ${FILTER_NONBAC})
  -h, --help                 Show this help message.
EOF
    exit 0
}

# --- Parse Arguments ---
ARGS=$(getopt -o d:t:h --long dbdir:,threads:,raw-reads-dir:,metadata-file:,temp-dir:,output-dir:,otutab-filename:,otus-filename:,primer-forward:,primer-reverse:,cutadapt-error-rate:,min-len:,max-len:,fastq-maxee:,fastq-qmax:,fastq_qmax:,fastq-max:,fastq_max:,min-unique-size:,id-cutoff:,sintax-db:,sintax-cutoff:,filter-nonbac:,help -n "$0" -- "$@")
if [ $? != 0 ]; then
    echo "Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"

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
        --primer-forward) PRIMER_FORWARD="$2"; shift 2 ;;
        --primer-reverse) PRIMER_REVERSE="$2"; shift 2 ;;
        --cutadapt-error-rate) CUTADAPT_ERROR_RATE="$2"; shift 2 ;;
        --min-len) MIN_LEN="$2"; shift 2 ;;
        --max-len) MAX_LEN="$2"; shift 2 ;;
        --fastq-maxee) FASTQ_MAXEE="$2"; shift 2 ;;
        --fastq-qmax|--fastq_qmax|--fastq-max|--fastq_max) FASTQ_QMAX="$2"; shift 2 ;;
        --min-unique-size) MIN_UNIQUE_SIZE="$2"; shift 2 ;;
        --id-cutoff) ID_CUTOFF="$2"; shift 2 ;;
        --sintax-db) SINTAX_DB="$2"; shift 2 ;;
        --sintax-cutoff) SINTAX_CUTOFF="$2"; shift 2 ;;
        --filter-nonbac) FILTER_NONBAC="$2"; shift 2 ;;
        -h|--help) usage ;;
        --) shift; break ;;
        *) echo "Internal error!" >&2; exit 1 ;;
    esac
done

: "${SINTAX_DB:=${DBDIR}/usearch/SILVA_modified.fasta}"
FILTER_NONBAC_R="${DBDIR}/script/otutab_filter_nonBac.R"
if [ "${FILTER_NONBAC}" = "TRUE" ] && [ ! -f "${FILTER_NONBAC_R}" ]; then
    echo "Error: non-bacteria filter script not found: ${FILTER_NONBAC_R}" >&2
    echo "  Use --dbdir /path/to/EasyMicrobiome  or  export EASY_MICROBIOME=..." >&2
    echo "  Or  --filter-nonbac FALSE" >&2
    exit 1
fi

ensure_uncompressed "${SINTAX_DB}" "SINTAX taxonomy database"
require_file "${METADATA_FILE}" "metadata file"

# --- Path Setup ---
MERGED_READS_DIR="${TEMP_DIR}"
CONCAT_MERGED_READS_FILE="${TEMP_DIR}/all.fq"
TRIMMED_FASTQ_FILE="${TEMP_DIR}/allfilter.fq"
FILTERED_READS_FILE="${TEMP_DIR}/filtered.fa"
DEREPLICATED_SEQS_FILE="${TEMP_DIR}/uniques.fa"
ASV_OTU_SEQS_FILE="${TEMP_DIR}/otus.fa"

RAW_OTU_TABLE_FILE="${OUTPUT_DIR}/raw/otutab.txt"
RAW_TAXONOMY_FILE="${OUTPUT_DIR}/raw/otus.sintax"
RAW_OTU_SEQS_FILE="${OUTPUT_DIR}/raw/otus.fa"

FILTERED_OTU_TABLE_FILE="${OUTPUT_DIR}/${OTUTAB_FILENAME}"
FILTERED_OTU_SEQS_FILE="${OUTPUT_DIR}/${OTUS_FILENAME}"
FILTERED_TAXONOMY_FILE="${OUTPUT_DIR}/otus.sintax"

mkdir -p "${TEMP_DIR}" "${OUTPUT_DIR}" "${OUTPUT_DIR}/raw"

for cmd in vsearch cutadapt Rscript csvtk gzip; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "Error: '${cmd}' not found. Please run in your EasyAmplicon conda environment." >&2
        exit 1
    fi
done

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

ensure_uncompressed() {
    local target="$1"
    local what="${2:-file}"
    if [ -f "$target" ]; then
        return 0
    fi
    if [ -f "${target}.gz" ]; then
        echo "Decompressing ${target}.gz -> ${target}"
        gzip -dc "${target}.gz" > "${target}"
        return 0
    fi
    echo "Error: missing ${what}: ${target} (or ${target}.gz)" >&2
    echo "  Download link (TODO): <PUT_DOWNLOAD_LINK_HERE>" >&2
    exit 1
}

echo "--- [feature_nanopore] Reads rename ---"
for i in $(tail -n+2 "${METADATA_FILE}" | cut -f1); do
    input_file=""
    for ext in fastq fq fastq.gz fq.gz; do
        if [ -f "${RAW_READS_DIR}/${i}.${ext}" ]; then
            input_file="${RAW_READS_DIR}/${i}.${ext}"
            break
        fi
    done

    if [ -z "${input_file}" ]; then
        echo "Error: cannot find reads for sample '${i}' under ${RAW_READS_DIR} (.fastq/.fq/.fastq.gz/.fq.gz)." >&2
        exit 1
    fi

    if [[ "${input_file}" == *.gz ]]; then
        decompressed_fastq="${TEMP_DIR}/${i}.decompressed.fastq"
        gzip -dc "${input_file}" > "${decompressed_fastq}"
        relabel_input="${decompressed_fastq}"
    else
        relabel_input="${input_file}"
    fi

    # vsearch doesn't support all usearch relabel behaviors, so we use --fastx_filter to relabel
    vsearch --fastx_filter "${relabel_input}" \
      --fastq_qmax "${FASTQ_QMAX}" \
      --fastqout "${MERGED_READS_DIR}/${i}.merged.fastq" \
      --relabel "${i}."
done

echo "--- [feature_nanopore] Integrate renamed reads ---"
cat "${MERGED_READS_DIR}"/*.merged.fastq > "${CONCAT_MERGED_READS_FILE}"
require_nonempty "${CONCAT_MERGED_READS_FILE}" "integrated merged FASTQ"
ls -lsh "${CONCAT_MERGED_READS_FILE}"

echo "--- [feature_nanopore] Primer trim and quality filter ---"
cutadapt -g "${PRIMER_FORWARD}...${PRIMER_REVERSE}" \
  --error-rate="${CUTADAPT_ERROR_RATE}" \
  --action=trim --rc \
  -j "${THREADS}" \
  --discard-untrimmed \
  -o "${TRIMMED_FASTQ_FILE}" \
  "${CONCAT_MERGED_READS_FILE}"
require_nonempty "${TRIMMED_FASTQ_FILE}" "cutadapt-trimmed FASTQ"

vsearch --fastx_filter "${TRIMMED_FASTQ_FILE}" \
  --fastq_minlen "${MIN_LEN}" \
  --fastq_maxlen "${MAX_LEN}" \
  --fastq_maxee "${FASTQ_MAXEE}" \
  --fastq_qmax "${FASTQ_QMAX}" \
  --fastaout "${FILTERED_READS_FILE}"
require_nonempty "${FILTERED_READS_FILE}" "filtered reads FASTA"

echo "--- [feature_nanopore] Dereplicate and cluster OTUs ---"
vsearch --derep_fulllength "${FILTERED_READS_FILE}" \
  --fasta_width 0 \
  --sizeout \
  --relabel Uni_ \
  --minuniquesize "${MIN_UNIQUE_SIZE}" \
  --threads "${THREADS}" \
  --output "${DEREPLICATED_SEQS_FILE}"
require_nonempty "${DEREPLICATED_SEQS_FILE}" "dereplicated sequences FASTA"

vsearch --cluster_unoise "${DEREPLICATED_SEQS_FILE}" \
  --minsize "${MIN_UNIQUE_SIZE}" \
  --centroids "${ASV_OTU_SEQS_FILE}"
require_nonempty "${ASV_OTU_SEQS_FILE}" "denoised centroids FASTA"

# Normalize centroid IDs so they match feature table IDs and R filtering logic.
# - Convert Uni_/Centroid_/Zotu_ to ASV_
# - Remove ";size=..." suffix to avoid mismatch with otutab rownames.
ASV_OTU_SEQS_FILE_NOSIZE="${TEMP_DIR}/otus_nosize.fa"
sed -E '
  s/^>Zotu(_)?/>ASV_/;
  s/^>Centroid_/>ASV_/;
  s/^>Uni_/>ASV_/;
  s/^>([^;[:space:]]+);size=[0-9]+.*/>\1/;
' "${ASV_OTU_SEQS_FILE}" > "${ASV_OTU_SEQS_FILE_NOSIZE}"
require_nonempty "${ASV_OTU_SEQS_FILE_NOSIZE}" "feature sequences FASTA (nosize)"

echo "--- [feature_nanopore] Reference chimera filtering ---"
vsearch --uchime_ref "${ASV_OTU_SEQS_FILE_NOSIZE}" \
  -db "${SINTAX_DB}" \
  --nonchimeras "${RAW_OTU_SEQS_FILE}" \
  --threads "${THREADS}"
sed -i 's/\r//g' "${RAW_OTU_SEQS_FILE}"
require_nonempty "${RAW_OTU_SEQS_FILE}" "chimera-filtered feature sequences FASTA"

echo "--- [feature_nanopore] Build raw feature table ---"
vsearch --usearch_global "${FILTERED_READS_FILE}" \
  --db "${RAW_OTU_SEQS_FILE}" \
  --id "${ID_CUTOFF}" \
  --threads "${THREADS}" \
  --otutabout "${RAW_OTU_TABLE_FILE}"
sed -i 's/\r//g' "${RAW_OTU_TABLE_FILE}"
require_nonempty "${RAW_OTU_TABLE_FILE}" "raw feature table"
csvtk -t stat "${RAW_OTU_TABLE_FILE}"

echo "--- [feature_nanopore] Taxonomy annotation ---"
vsearch --sintax "${RAW_OTU_SEQS_FILE}" \
  --db "${SINTAX_DB}" \
  --sintax_cutoff "${SINTAX_CUTOFF}" \
  --tabbedout "${RAW_TAXONOMY_FILE}"
sed -i 's/\r//g' "${RAW_TAXONOMY_FILE}"
require_nonempty "${RAW_TAXONOMY_FILE}" "raw taxonomy table"

echo "--- [feature_nanopore] Filter contaminants (optional) ---"
if [ "${FILTER_NONBAC}" = "TRUE" ]; then
    Rscript "${FILTER_NONBAC_R}" \
      --input "${RAW_OTU_TABLE_FILE}" \
      --taxonomy "${RAW_TAXONOMY_FILE}" \
      --output "${FILTERED_OTU_TABLE_FILE}" \
      --stat "${OUTPUT_DIR}/raw/otutab_nonBac.stat" \
      --discard "${OUTPUT_DIR}/raw/otus.sintax.discard"
    require_nonempty "${FILTERED_OTU_TABLE_FILE}" "filtered feature table"

    cut -f 1 "${FILTERED_OTU_TABLE_FILE}" | tail -n+2 > "${OUTPUT_DIR}/otutab.id"
    require_nonempty "${OUTPUT_DIR}/otutab.id" "filtered feature ID list"
    vsearch --fastx_getseqs "${RAW_OTU_SEQS_FILE}" \
      --labels "${OUTPUT_DIR}/otutab.id" \
      --fastaout "${FILTERED_OTU_SEQS_FILE}"
    require_nonempty "${FILTERED_OTU_SEQS_FILE}" "filtered feature sequences FASTA"
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}' \
      "${RAW_TAXONOMY_FILE}" "${OUTPUT_DIR}/otutab.id" \
      > "${FILTERED_TAXONOMY_FILE}"
    require_nonempty "${FILTERED_TAXONOMY_FILE}" "filtered taxonomy table"
else
    cp -f "${RAW_OTU_TABLE_FILE}" "${FILTERED_OTU_TABLE_FILE}"
    cp -f "${RAW_OTU_SEQS_FILE}" "${FILTERED_OTU_SEQS_FILE}"
    cp -f "${RAW_TAXONOMY_FILE}" "${FILTERED_TAXONOMY_FILE}"
    require_nonempty "${FILTERED_OTU_TABLE_FILE}" "feature table"
    require_nonempty "${FILTERED_OTU_SEQS_FILE}" "feature sequences FASTA"
    require_nonempty "${FILTERED_TAXONOMY_FILE}" "taxonomy table"
fi

echo "--- [feature_nanopore] Done ---"