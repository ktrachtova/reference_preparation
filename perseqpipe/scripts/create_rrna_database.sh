#!/bin/bash
# Script to create a custom rRNA database

# Set up project directory
PROJECT_DIR=$(dirname "$(pwd)")
echo "Project directory: $PROJECT_DIR"

CREATE_FASTA_MMSEQS2=$(pwd)/utils/create_fasta_mmseqs2.py

DB1="" # RNACentral
DB2="" # NCBI

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --DB1) DB1="$2"; shift ;;
        --DB2) DB2="$2"; shift ;;
        --help)
            echo "Usage: $0 --DB1 <file> --DB2 <file>"
            exit 0 ;;
        *)
            echo "Unknown parameter passed: $1"
            exit 1 ;;
    esac
    shift
done

# Check if both DB1 and DB2 are supplied
if [[ -z "$DB1" || -z "$DB2" ]]; then
    echo "Error: Both --DB1 and --DB2 must be provided."
    echo "Usage: $0 --DB1 <file> --DB2 <file>"
    exit 1
fi

# Decompress if gzipped
if [[ "$DB1" == *.gz ]]; then
    gzip -d -c "$DB1" > "${DB1%.gz}"
    DB1="${DB1%.gz}"
fi

if [[ "$DB2" == *.gz ]]; then
    gzip -d -c "$DB2" > "${DB2%.gz}"
    DB2="${DB2%.gz}"
fi

echo "Using RNACentral rRNA FASTA file: $DB1"
echo "Using NCBI rRNA FASTA file: $DB2"

# Output directory for results
OUTPUT_DIR="${PROJECT_DIR}/reference_files/rRNA/$(date +'%Y_%m_%d')"
mkdir -p "$OUTPUT_DIR"

# Create temporary folder for intermediate files
TMP_DIR="$(pwd)/tmp"
mkdir -p $TMP_DIR

echo ""
echo "--------------------------"
echo "Cleaning rRNA sequences..."

# Convert FASTA to one-line format and replace Us with Ts
process_fasta() {
    local input=$1
    local output=$2
    awk '/^>/ {printf("\n%s\n",$0);next;} {printf("%s",$0);} END {printf("\n");}' < "$input" | sed '1d' | grep "\S" |
        awk '{ if ($0 !~ /^>/) gsub(/U/, "T"); print }' > "$output"
}

process_fasta "$DB1" ${TMP_DIR}/"_db1_cleaned.fa"
process_fasta "$DB2" ${TMP_DIR}/"_db2_cleaned.fa"

cat ${TMP_DIR}/_db1_cleaned.fa ${TMP_DIR}/_db2_cleaned.fa > ${TMP_DIR}/tmp.fa

#docker run --rm \
#    -v "$(pwd)/tmp.fa:/data/tmp.fa" \
#    -v "$OUTPUT_DIR:/output" \
#    ktrachtok/reference_preparation:latest \
#    /cd-hit-v4.8.1-2019-0228/cd-hit-est -i /data/tmp.fa -c 1 -s 1 -aL 1 -aS 1 -d 0 -p 1 -g 1 -o /output/rRNA_db_custom.fa

echo ""
echo "-------------------------------"
echo "Clustering with MMseqs2..."

docker run --rm \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:latest \
    bash -c "
        /MMseqs2/build/bin/mmseqs createdb /data/tmp.fa /data/inputDB && \
        /MMseqs2/build/bin/mmseqs clusthash /data/inputDB /data/resultDB --min-seq-id 1.0 && \
        /MMseqs2/build/bin/mmseqs clust /data/inputDB /data/resultDB /data/clusterDB && \
        /MMseqs2/build/bin/mmseqs createtsv /data/inputDB /data/inputDB /data/clusterDB /data/cluster_result.tsv
    "

echo ""
echo "-------------------------------"
echo "Generating clustered FASTA file..."

docker run --rm \
    -v "${TMP_DIR}:/data" \
    -v ${CREATE_FASTA_MMSEQS2}:/scripts/create_fasta_mmseqs2.py \
    ktrachtok/reference_preparation:latest \
    python3 /scripts/create_fasta_mmseqs2.py \
        --fasta /data/tmp.fa \
        --mmseqs2_tsv /data/cluster_result.tsv \
        --output /data/rRNA_db_custom.fa


# Cleanup
# rm tmp.fa _*

echo ""
echo "-------------------------------"
echo "Calculating final statistics..."

MIN_LENGTH=$(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' "$TMP_DIR/rRNA_db_custom.fa" | sort -n | head -n 1)
MAX_LENGTH=$(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' "$TMP_DIR/rRNA_db_custom.fa" | sort -n | tail -n 1)
AVG_LENGTH=$(awk '/^>/ {if (seqlen) {sum += seqlen; count++} seqlen=0; next} {seqlen += length($0)} END {print sum/count}' "$TMP_DIR/rRNA_db_custom.fa")
NUM_SEQ=$(grep -c "^>" "$TMP_DIR/rRNA_db_custom.fa")

echo "Number of rRNA sequences: $NUM_SEQ"
echo "Minimum length: $MIN_LENGTH"
echo "Maximum length: $MAX_LENGTH"
echo "Average length: $AVG_LENGTH"

cp ${TMP_DIR}/rRNA_db_custom.fa ${OUTPUT_DIR}/rRNA_db_custom.fa

# Cleaning
rm -rf $TMP_DIR