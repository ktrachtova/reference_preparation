#!/bin/bash
set -euo pipefail

source ../config.mk

# Set up project directory
PROJECT_DIR=$(dirname "$(pwd)")
mkdir -p "$PROJECT_DIR/databases"

echo "Project directory: $PROJECT_DIR"
echo "Docker image version: $REFERENCE_PREPARATION_DOCKER_VERSION"

# -----------------------------
# Script: create_custom_tRNA_db.sh
# Purpose: Build a custom tRNA database
# -----------------------------

# Inputs
DB1=""           # GENCODE FASTA file
DB2=""           # GtRNAdb BED file
DB2_FASTA1=""    # filtered-tRNAs.fa
DB2_FASTA2=""    # mature-tRNAs.fa
GENOME=""        # Genome FASTA file (optional)

# PSL to BED12 R conversion scripts
PSL2BED=$(pwd)/utils/psl2bed.r
BED2GTF=$(pwd)/utils/bed2gtf.py

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --DB1) DB1="$2"; shift ;;
        --DB2) DB2="$2"; shift ;;
        --DB2_FASTA1) DB2_FASTA1="$2"; shift ;;
        --DB2_FASTA2) DB2_FASTA2="$2"; shift ;;
        --GENOME) GENOME="$2"; shift ;;
        --help)
            echo "Usage: $0 --DB1 <gencode.fa> --DB2 <trna.bed> --DB2_FASTA1 <filtered.fa> --DB2_FASTA2 <mature.fa> [--GENOME <genome.fa>]"
            exit 0 ;;
        *) echo "❌ Unknown parameter: $1" >&2; exit 1 ;;
    esac
    shift
done

# Validate inputs
if [[ -z "$DB1" || -z "$DB2" || -z "$DB2_FASTA1" || -z "$DB2_FASTA2" ]]; then
    echo "❌ All input files (except --GENOME) are required." >&2
    echo "Use --help for usage." >&2
    exit 1
fi

# Decompress GENCODE if needed
if [[ "$DB1" == *.gz ]]; then
    gzip -d -c "$DB1" > "${DB1%.gz}"
    DB1="${DB1%.gz}"
fi

# Decompress DB2_FASTA files if needed
for var in DB2_FASTA1 DB2_FASTA2 DB2; do
    val="${!var}"
    if [[ "$val" == *.gz ]]; then
        gzip -d -c "$val" > "${val%.gz}"
        eval "$var='${val%.gz}'"
    fi
    echo "Using $(echo "$var" | tr '[:lower:]' '[:upper:]'): ${!var}"
done

# Genome file setup
if [[ -z "$GENOME" ]]; then
    echo "🧬 Genome FASTA not provided. Downloading..."
    wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz -O "$PROJECT_DIR/databases/genome.fa.gz"
    gunzip -c "$PROJECT_DIR/databases/genome.fa.gz" > "$PROJECT_DIR/databases/genome.fa"
    GENOME="$PROJECT_DIR/databases/genome.fa"
else
    if [[ "$GENOME" == *.gz ]]; then
        gunzip -c "$GENOME" > "${GENOME%.gz}"
        GENOME="${GENOME%.gz}"
    fi
fi


# Output directory
OUTPUT_DIR="$PROJECT_DIR/reference_files/tRNA/$(date +'%Y_%m_%d')"
mkdir -p "$OUTPUT_DIR"

# Create temporary folder for intermediate files
TMP_DIR="$(pwd)/tmp"
mkdir -p $TMP_DIR

echo ""
echo "--------------------------"
echo "Cleaning tRNA sequences..."

# FASTA multiline to one line for Gencode
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' "$DB1" > ${TMP_DIR}/gencode.mt_tRNA.oneline.fa

# Extract only Mt-tRNA sequences from GENCODE
grep -A 1 "Mt_tRNA" ${TMP_DIR}/gencode.mt_tRNA.oneline.fa | grep -v "^--" > $TMP_DIR/mt_tRNA_seq.fa

# The Mt-tRNA sequences contain "|" in their headers, we need to replace the pipe with "_"
# as the "|" is a separator for IDs from multiple databases for the same sequence in the sncRNA GTF file;
# see https://github.com/ktrachtova/perseqpipe/blob/main/docs/annotation_preparation.md#sncrna-gtf-file-format
sed '/^>/{ s/|$//; s/|/_/g; }' "$TMP_DIR/mt_tRNA_seq.fa" > "$TMP_DIR/mt_tRNA_db.fa"

# FASTA multiline to one line for GtRNAdb
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' "$DB2_FASTA1" > ${TMP_DIR}/tmp1_
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' "$DB2_FASTA2" > ${TMP_DIR}/tmp2_

# Replace Us with Ts
awk 'BEGIN { OFS="\n" } \
    /^>/ { print; next } \
    { gsub(/U/, "T"); gsub(/u/, "t"); print }' "${TMP_DIR}"/mt_tRNA_db.fa > "${TMP_DIR}"/mt_tRNA_db.fa.tmp
mv "${TMP_DIR}"/mt_tRNA_db.fa.tmp "${TMP_DIR}"/mt_tRNA_db.fa

awk 'BEGIN { OFS="\n" } \
    /^>/ { print; next } \
    { gsub(/U/, "T"); gsub(/u/, "t"); print }' ${TMP_DIR}/tmp1_ > ${TMP_DIR}/tmp1_.tmp
mv ${TMP_DIR}/tmp1_.tmp ${TMP_DIR}/tmp1_

awk 'BEGIN { OFS="\n" } \
    /^>/ { print; next } \
    { gsub(/U/, "T"); gsub(/u/, "t"); print }' ${TMP_DIR}/tmp2_ > ${TMP_DIR}/tmp2_.tmp
mv ${TMP_DIR}/tmp2_.tmp ${TMP_DIR}/tmp2_

# Map Mt-tRNA to genome wtih BLAT
echo ""
echo "------------------------------------"
echo "Aligning tRNA sequences to genome..."

docker run --rm \
    -v "${TMP_DIR}:/data" \
    -v $GENOME:/data/genome.fa \
    ktrachtok/reference_preparation:x86_64-"${REFERENCE_PREPARATION_DOCKER_VERSION}" \
    ./blat /data/genome.fa /data/mt_tRNA_db.fa \
    -t=dna -q=rna -maxIntron=10000 -stepSize=5 -repMatch=1000 \
    -minScore=20 -minIdentity=100 -noTrimA -out=psl \
    /data/mt_tRNA_db_genomeMap.psl

# Convert PSL format to BED12 format
echo ""
echo "---------------------"
echo "Converting PSL to BED"

docker run --rm \
    -v ${PSL2BED}:/scripts/psl2bed.r \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:x86_64-"${REFERENCE_PREPARATION_DOCKER_VERSION}" \
    Rscript /scripts/psl2bed.r /data/mt_tRNA_db_genomeMap.psl /data/mt_tRNA_db_genomeMap.bed

# Merge GtfRNAdb BED12 file with tRNA with created BED12 file with Mt-tRNA
cat ${TMP_DIR}/mt_tRNA_db_genomeMap.bed $DB2 > ${TMP_DIR}/tRNA_db_custom_genomeMap.bed

# Convert BED12 to GTF
echo ""
echo "----------------------"
echo "Converting BED to GTF"

docker run --rm \
    -v ${BED2GTF}:/scripts/bed2gtf.py \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:x86_64-"${REFERENCE_PREPARATION_DOCKER_VERSION}" \
    python3 /scripts/bed2gtf.py -i /data/tRNA_db_custom_genomeMap.bed -o /data/tRNA_db_custom_genomeMap.gtf --gene_feature --gene_biotype tRNA --source GtRNAdb

# Length distribution of final tRNA database
echo ""
echo "--------------------------"
echo "Making final statistics..."

echo "tRNA_length tRNA_number" > ${OUTPUT_DIR}/trna_database_lenDist.csv
cat ${TMP_DIR}/mt_tRNA_db.fa ${TMP_DIR}/tmp1_ ${TMP_DIR}/tmp2_ > ${TMP_DIR}/tRNA_db_custom.fa

cat ${TMP_DIR}/tRNA_db_custom.fa | awk 'NR%2 == 0 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' >> ${OUTPUT_DIR}/trna_database_lenDist.csv

# Calculate and print some statistics about created Mt-tRNA database
# Calculate minimum, maximum, and average lengths
read -r MIN_LENGTH MAX_LENGTH AVG_LENGTH NUM_SEQ < <(
awk '
  function flush() {
    if (seqlen>0) {
      if (count==0 || seqlen<min) min=seqlen
      if (count==0 || seqlen>max) max=seqlen
      sum+=seqlen
      count++
      seqlen=0
    }
  }
  /^>/ { flush(); next }
  { seqlen += length($0) }
  END {
    flush()
    avg = (count? sum/count : 0)
    printf "%d %d %.6f %d\n", min, max, avg, count
  }
' "$TMP_DIR/tRNA_db_custom.fa"
)

# Print the results
echo "Number of tRNA + Mt-tRNA sequences: $NUM_SEQ"
echo "Minimum length: $MIN_LENGTH"
echo "Maximum length: $MAX_LENGTH"
echo "Average length: $AVG_LENGTH"

cp ${TMP_DIR}/tRNA_db_custom.fa ${OUTPUT_DIR}/tRNA_db_custom.fa
cp ${TMP_DIR}/tRNA_db_custom_genomeMap.gtf ${OUTPUT_DIR}/tRNA_db_custom_genomeMap.gtf
cp ${TMP_DIR}/tRNA_db_custom_genomeMap.bed ${OUTPUT_DIR}/tRNA_db_custom_genomeMap.bed
cp ${TMP_DIR}/mt_tRNA_db_genomeMap.psl ${OUTPUT_DIR}/mt_tRNA_db_genomeMap.psl
cp ${TMP_DIR}/mt_tRNA_db_genomeMap.bed ${OUTPUT_DIR}/mt_tRNA_db_genomeMap.bed

# Cleaning
rm -rf $TMP_DIR
