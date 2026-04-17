#!/bin/bash
set -e
# Module: feature_pacbio - PacBio feature table construction (stop before normalization)

# --- Default Parameters ---
# EasyMicrobiome root (contains script/, usearch/, etc.). Override per machine:
#   export EASY_MICROBIOME=/path/to/EasyMicrobiome
#   or: bash easyamplicon.sh feature_pacbio --dbdir /path/to/EasyMicrobiome ...
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
FASTQ_MAXEE=2
FASTQ_QMAX=93
MIN_UNIQUE_SIZE=2
ID_CUTOFF=0.97
SINTAX_DB=""
SINTAX_CUTOFF=0.1
FILTER_NONBAC="TRUE"

usage() {
    cat <<EOF
Module: feature_pacbio - PacBio feature table construction (stops before normalization)

Function: Relabel reads, cutadapt primer trim, length/quality filter, dereplicate,
          UNOISE/ASV (usearch or vsearch fallback), feature table, SINTAX, optional non-bac filter.

Usage: bash easyamplicon.sh feature_pacbio [OPTIONS]

Options:
  -d, --dbdir DIR            EasyMicrobiome installation root (script/, usearch/).
                             Default: \$EASY_MICROBIOME if set, else /d/EasyMicrobiome.
                             (Current default: ${DBDIR})
  -t, --threads INT          Threads for parallel steps (cutadapt, vsearch).
                             (Default: ${THREADS})

  --raw-reads-dir DIR        Directory with per-sample reads: <SampleID>.fastq|.fq|.fastq.gz|.fq.gz
                             (Default: ${RAW_READS_DIR})
  --metadata-file FILE       Metadata TSV; column 1 = SampleID (matches filenames).
                             (Default: ${METADATA_FILE})
  --temp-dir DIR             Temporary/intermediate files.
                             (Default: ${TEMP_DIR})
  --output-dir DIR           Results root (creates raw/ and final tables).
                             (Default: ${OUTPUT_DIR})
  --otutab-filename FILE     Filtered feature table filename under output-dir.
                             (Default: ${OTUTAB_FILENAME})
  --otus-filename FILE       Filtered representative sequences filename under output-dir.
                             (Default: ${OTUS_FILENAME})

  --primer-forward STR       Forward primer (cutadapt -g F...R).
                             (Default: ${PRIMER_FORWARD})
  --primer-reverse STR       Reverse primer.
                             (Default: ${PRIMER_REVERSE})
  --cutadapt-error-rate F    Cutadapt --error-rate.
                             (Default: ${CUTADAPT_ERROR_RATE})
  --min-len INT              vsearch --fastq_minlen after cutadapt.
                             (Default: ${MIN_LEN})
  --max-len INT              vsearch --fastq_maxlen.
                             (Default: ${MAX_LEN})
  --fastq-maxee FLOAT        vsearch --fastq_maxee on trimmed reads.
                             (Default: ${FASTQ_MAXEE})
  --fastq-qmax, --fastq_qmax, --fastq-max, --fastq_max INT
                             FASTQ quality ceiling (PacBio/ONT often 93).
                             (Default: ${FASTQ_QMAX})
  --min-unique-size INT      Dereplication / UNOISE minsize.
                             (Default: ${MIN_UNIQUE_SIZE})
  --id-cutoff FLOAT          vsearch --usearch_global --id.
                             (Default: ${ID_CUTOFF})
  --sintax-db FILE           SINTAX reference FASTA.
                             (Default: ${SINTAX_DB})
  --sintax-cutoff FLOAT      SINTAX confidence cutoff.
                             (Default: ${SINTAX_CUTOFF})
  --filter-nonbac BOOL       Run otutab_filter_nonBac.R (TRUE/FALSE).
                             (Default: ${FILTER_NONBAC})

  -h, --help                 Show this help and exit.

Note: usearch is optional; if not runnable, UNOISE uses vsearch --cluster_unoise.
EOF
    exit 0
}

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

: "${SINTAX_DB:=${DBDIR}/usearch/sintax_defalut_emu_database.fasta}"

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

ensure_uncompressed() {
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

FILTER_NONBAC_R="${DBDIR}/script/otutab_filter_nonBac.R"
if [ "${FILTER_NONBAC}" = "TRUE" ] && [ ! -f "${FILTER_NONBAC_R}" ]; then
    echo "Error: non-bacteria filter script not found: ${FILTER_NONBAC_R}" >&2
    echo "  Point to your EasyMicrobiome root, e.g.  --dbdir /mnt/d/EasyMicrobiome" >&2
    echo "  Or set once:  export EASY_MICROBIOME=/path/to/EasyMicrobiome" >&2
    echo "  Or skip this step:  --filter-nonbac FALSE" >&2
    exit 1
fi

ensure_uncompressed "${SINTAX_DB}" "SINTAX taxonomy database"
require_file "${METADATA_FILE}" "metadata file"

for cmd in vsearch cutadapt Rscript csvtk; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "Error: '${cmd}' not found. Please run in your EasyAmplicon conda environment." >&2
        exit 1
    fi
done

HAS_USEARCH="TRUE"
if ! command -v usearch >/dev/null 2>&1; then
    HAS_USEARCH="FALSE"
elif ! usearch -version >/dev/null 2>&1; then
    HAS_USEARCH="FALSE"
fi

MERGED_READS_DIR="${TEMP_DIR}"
CONCAT_MERGED_READS_FILE="${TEMP_DIR}/all.fq"
TRIMMED_FASTQ_FILE="${TEMP_DIR}/allfilter.fq"
FILTERED_READS_FILE="${TEMP_DIR}/filtered.fa"
DEREPLICATED_SEQS_FILE="${TEMP_DIR}/uniques.fa"
ZOTUS_FILE="${TEMP_DIR}/zotus.fa"
ASV_OTU_SEQS_FILE="${TEMP_DIR}/otus.fa"
RAW_OTU_TABLE_FILE="${OUTPUT_DIR}/raw/otutab.txt"
RAW_TAXONOMY_FILE="${OUTPUT_DIR}/raw/otus.sintax"
RAW_OTU_SEQS_FILE="${OUTPUT_DIR}/raw/otus.fa"
FILTERED_OTU_TABLE_FILE="${OUTPUT_DIR}/${OTUTAB_FILENAME}"
FILTERED_OTU_SEQS_FILE="${OUTPUT_DIR}/${OTUS_FILENAME}"
FILTERED_TAXONOMY_FILE="${OUTPUT_DIR}/otus.sintax"

mkdir -p "${TEMP_DIR}" "${OUTPUT_DIR}" "${OUTPUT_DIR}/raw"

echo "--- [feature_pacbio] Reads rename ---"
for i in $(tail -n+2 "${METADATA_FILE}" | cut -f1); do
    input_file=""
    for ext in fastq fq fastq.gz fq.gz; do
        if [ -f "${RAW_READS_DIR}/${i}.${ext}" ]; then
            input_file="${RAW_READS_DIR}/${i}.${ext}"
            break
        fi
    done
    if [ -z "${input_file}" ]; then
        echo "Error: cannot find reads for sample '${i}' under ${RAW_READS_DIR}." >&2
        exit 1
    fi
    if [[ "${input_file}" == *.gz ]]; then
        decompressed_fastq="${TEMP_DIR}/${i}.decompressed.fastq"
        gzip -dc "${input_file}" > "${decompressed_fastq}"
        relabel_input="${decompressed_fastq}"
    else
        relabel_input="${input_file}"
    fi

    vsearch --fastx_filter "${relabel_input}" \
      --fastq_qmax "${FASTQ_QMAX}" \
      --fastqout "${MERGED_READS_DIR}/${i}.merged.fastq" \
      --relabel "${i}."
done

echo "--- [feature_pacbio] Integrate renamed reads ---"
cat "${MERGED_READS_DIR}"/*.merged.fastq > "${CONCAT_MERGED_READS_FILE}"
require_nonempty "${CONCAT_MERGED_READS_FILE}" "integrated merged FASTQ"
ls -lsh "${CONCAT_MERGED_READS_FILE}"

echo "--- [feature_pacbio] Primer trim and quality filter ---"
require_file "${CONCAT_MERGED_READS_FILE}" "integrated FASTQ"
cutadapt -g "${PRIMER_FORWARD}...${PRIMER_REVERSE}" --error-rate="${CUTADAPT_ERROR_RATE}" --action=trim --rc -j "${THREADS}" --discard-untrimmed -o "${TRIMMED_FASTQ_FILE}" "${CONCAT_MERGED_READS_FILE}"
require_nonempty "${TRIMMED_FASTQ_FILE}" "cutadapt-trimmed FASTQ"

n_after_cutadapt=$(awk 'NR%4==1{c++} END{print c+0}' "${TRIMMED_FASTQ_FILE}")
if [ "${n_after_cutadapt}" -lt 500 ]; then
    echo "Warning: only ${n_after_cutadapt} reads passed cutadapt (linked primers)." >&2
    echo "  If this is unexpected, your --primer-forward/--primer-reverse may not match the library." >&2
    echo "  Tutorial PacBio often uses e.g. AGAGTTTGATCCTGGCTCAG...AAGTCSTAACAAGGTADCCSTA (not 1492R)." >&2
fi

vsearch --fastx_filter "${TRIMMED_FASTQ_FILE}" --fastq_minlen "${MIN_LEN}" --fastq_maxlen "${MAX_LEN}" --fastq_maxee "${FASTQ_MAXEE}" --fastq_qmax "${FASTQ_QMAX}" --fastaout "${FILTERED_READS_FILE}"
require_nonempty "${FILTERED_READS_FILE}" "filtered reads FASTA"

n_fasta=$(grep -c '^>' "${FILTERED_READS_FILE}" 2>/dev/null || echo 0)
if [ "${n_fasta}" -eq 0 ]; then
    echo "Error: zero sequences in ${FILTERED_READS_FILE} after length/quality filter." >&2
    echo "  Often caused by cutadapt leaving almost no reads (check primers) or MIN_LEN/MAX_LEN too strict." >&2
    exit 1
fi

echo "--- [feature_pacbio] Dereplicate and denoise ---"
vsearch --derep_fulllength "${FILTERED_READS_FILE}" --fasta_width 0 --sizeout --relabel Uni_ --minuniquesize "${MIN_UNIQUE_SIZE}" --output "${DEREPLICATED_SEQS_FILE}"
require_nonempty "${DEREPLICATED_SEQS_FILE}" "dereplicated sequences FASTA"

n_uniq=$(grep -c '^>' "${DEREPLICATED_SEQS_FILE}" 2>/dev/null || echo 0)
if [ "${n_uniq}" -eq 0 ]; then
    echo "Error: dereplication wrote zero sequences to ${DEREPLICATED_SEQS_FILE}." >&2
    echo "  With only ${n_fasta} input sequences, --min-unique-size ${MIN_UNIQUE_SIZE} may remove all singletons." >&2
    echo "  Try: --min-unique-size 1   and fix cutadapt/ primers so most reads survive trimming." >&2
    exit 1
fi

if [ "${HAS_USEARCH}" = "TRUE" ]; then
    usearch -unoise3 "${DEREPLICATED_SEQS_FILE}" -minsize "${MIN_UNIQUE_SIZE}" -zotus "${ZOTUS_FILE}"
else
    echo "Warning: usearch unavailable, fallback to vsearch --cluster_unoise." >&2
    vsearch --cluster_unoise "${DEREPLICATED_SEQS_FILE}" \
      --minsize "${MIN_UNIQUE_SIZE}" \
      --centroids "${ZOTUS_FILE}"
fi
require_nonempty "${ZOTUS_FILE}" "denoised centroids FASTA"
sed -E '
  s/^>Zotu(_)?/>ASV_/;
  s/^>Centroid_/>ASV_/;
  s/^>Uni_/>ASV_/;
' "${ZOTUS_FILE}" > "${ASV_OTU_SEQS_FILE}"
require_nonempty "${ASV_OTU_SEQS_FILE}" "ASV representative sequences FASTA"
ASV_OTU_SEQS_FILE_NOSIZE="${TEMP_DIR}/otus_nosize.fa"
# Strip ";size=..." to make OTU IDs match the feature table (which uses ASV_1, not ASV_1;size=...).
sed -E 's/^>([^;[:space:]]+);size=[0-9]+.*/>\1/' "${ASV_OTU_SEQS_FILE}" > "${ASV_OTU_SEQS_FILE_NOSIZE}"
require_nonempty "${ASV_OTU_SEQS_FILE_NOSIZE}" "ASV representative sequences FASTA (nosize)"
cp -f "${ASV_OTU_SEQS_FILE_NOSIZE}" "${RAW_OTU_SEQS_FILE}"
require_nonempty "${RAW_OTU_SEQS_FILE}" "raw feature sequences FASTA"

echo "--- [feature_pacbio] Build raw feature table ---"
vsearch --usearch_global "${FILTERED_READS_FILE}" --db "${RAW_OTU_SEQS_FILE}" --id "${ID_CUTOFF}" --threads "${THREADS}" --otutabout "${RAW_OTU_TABLE_FILE}"
require_nonempty "${RAW_OTU_TABLE_FILE}" "raw feature table"
sed -i 's/\r//g' "${RAW_OTU_TABLE_FILE}"
csvtk -t stat "${RAW_OTU_TABLE_FILE}"

echo "--- [feature_pacbio] Taxonomy annotation ---"
vsearch --sintax "${RAW_OTU_SEQS_FILE}" --db "${SINTAX_DB}" --sintax_cutoff "${SINTAX_CUTOFF}" --tabbedout "${RAW_TAXONOMY_FILE}"
sed -i 's/\r//g' "${RAW_TAXONOMY_FILE}"
require_nonempty "${RAW_TAXONOMY_FILE}" "raw taxonomy table"

echo "--- [feature_pacbio] Filter contaminants (optional) ---"
if [ "${FILTER_NONBAC}" = "TRUE" ]; then
    Rscript "${FILTER_NONBAC_R}" --input "${RAW_OTU_TABLE_FILE}" --taxonomy "${RAW_TAXONOMY_FILE}" --output "${FILTERED_OTU_TABLE_FILE}" --stat "${OUTPUT_DIR}/raw/otutab_nonBac.stat" --discard "${OUTPUT_DIR}/raw/otus.sintax.discard"
    require_nonempty "${FILTERED_OTU_TABLE_FILE}" "filtered feature table"
    cut -f 1 "${FILTERED_OTU_TABLE_FILE}" | tail -n+2 > "${OUTPUT_DIR}/otutab.id"
    require_nonempty "${OUTPUT_DIR}/otutab.id" "filtered feature ID list"
    vsearch --fastx_getseqs "${RAW_OTU_SEQS_FILE}" --labels "${OUTPUT_DIR}/otutab.id" --fastaout "${FILTERED_OTU_SEQS_FILE}"
    require_nonempty "${FILTERED_OTU_SEQS_FILE}" "filtered feature sequences FASTA"
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}' "${RAW_TAXONOMY_FILE}" "${OUTPUT_DIR}/otutab.id" > "${FILTERED_TAXONOMY_FILE}"
    require_nonempty "${FILTERED_TAXONOMY_FILE}" "filtered taxonomy table"
else
    cp -f "${RAW_OTU_TABLE_FILE}" "${FILTERED_OTU_TABLE_FILE}"
    cp -f "${RAW_OTU_SEQS_FILE}" "${FILTERED_OTU_SEQS_FILE}"
    cp -f "${RAW_TAXONOMY_FILE}" "${FILTERED_TAXONOMY_FILE}"
    require_nonempty "${FILTERED_OTU_TABLE_FILE}" "feature table"
    require_nonempty "${FILTERED_OTU_SEQS_FILE}" "feature sequences FASTA"
    require_nonempty "${FILTERED_TAXONOMY_FILE}" "taxonomy table"
fi

echo "--- [feature_pacbio] Done ---"
