#!/bin/bash
# Script to create a custom piRNA database
set -euo pipefail

source ../config.mk

# Set up project directory
PROJECT_DIR=$(dirname $(pwd))
echo "Project directory: $PROJECT_DIR"
mkdir -p $PROJECT_DIR # folder for downloaded reference files
echo "Docker image version: $REFERENCE_PREPARATION_DOCKER_VERSION"

# Inputs from command line
DB1="" # RNACentral FASTA
DB2="" # piRBase FASTA
DB3="" # piRNAdb FASTA
DB4="" # NCBI FASTA
GENOME="" # Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --DB1) DB1="$2"; shift ;;   # Capture value for --db1
        --DB2) DB2="$2"; shift ;;   # Capture value for --db2
        --DB3) DB3="$2"; shift ;;   # Capture value for --db3
        --DB4) DB4="$2"; shift ;;   # Capture value for --db4
        --GENOME) GENOME="$2"; shift ;;   # Capture value for --GENOME
        --help)  # Add a help flag
            echo "Usage: $0 --DB1 <file> --DB2 <file> --DB3 <file> --DB4 <file> --GENOME <file>"
            exit 0
            ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Validate required inputs
if [[ -z "$DB1" || -z "$DB2" || -z "$DB3" || -z "$DB4" ]]; then
    echo "Error: All 4 database FASTA files (--DB1, --DB2, --DB3, --DB4) must be provided."
    exit 1
fi

# Decompress db1–db4 if they are gzipped
for var in DB1 DB2 DB3 DB4; do
    val="${!var}"
    if [[ "$val" == *.gz ]]; then
        echo "Decompressing $val..."
        gzip -d -c "$val" > "${val%.gz}"
        eval "$var='${val%.gz}'"
    fi
    echo "Using $(echo "$var" | tr '[:lower:]' '[:upper:]'): ${!var}"
done

# Handle genome input
if [ -z "$GENOME" ]; then
    echo "🧬 Genome FASTA not provided. Downloading GRCh38..."
    if [ ! -f "$PROJECT_DIR/databases/GRCh38.primary_assembly.genome.fa" ]; then
        wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
        gzip -d GRCh38.primary_assembly.genome.fa.gz
        mv GRCh38.primary_assembly.genome.fa "$PROJECT_DIR/databases"
    fi
    GENOME="$PROJECT_DIR/databases/GRCh38.primary_assembly.genome.fa"
else
    if [[ "$GENOME" == *.gz ]]; then
        echo "Decompressing genome FASTA..."
        gunzip -c "$GENOME" > "${GENOME%.gz}"
        GENOME="${GENOME%.gz}"
    fi
    echo "Using genome FASTA: $GENOME"
fi

# R script to create sequence logos and PSL to BED12 conversion
PSL2BED=$(pwd)/utils/psl2bed.r
BED2GTF=$(pwd)/utils/bed2gtf.py
CREATE_FASTA_MMSEQS2=$(pwd)/utils/create_fasta_mmseqs2.py
# SEQUENCE_LOGO_R=$(pwd)/utils/sequence_logos.r

# Create temporary folder for intermediate files
TMP_DIR="$(pwd)/tmp"
mkdir -p $TMP_DIR

# Output directory
OUTPUT_DIR=${PROJECT_DIR}/reference_files/piRNA/$(date +'%Y_%m_%d')
echo "${OUTPUT_DIR}"
mkdir -p $OUTPUT_DIR

echo "Merging databases, converting to one-line FASTA..."
# Merge all piRNA databases and shorten name
cat $DB1 $DB2 $DB3 $DB3 | cut -d ' ' -f1 - > $TMP_DIR/_tmp

# In case FASTA is multiline, change to oneline
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' $TMP_DIR/_tmp > $TMP_DIR/tmp.oneline

# Remove sequences with Ns
echo ""
echo "-----------------------------"
echo "Removing sequences with Ns..."
docker run --rm \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:x86_64-"${REFERENCE_PREPARATION_DOCKER_VERSION}" \
    cutadapt --max-n 0 -o /data/tmp.filtered.fa /data/tmp.oneline

# Remove redundant sequences
echo ""
echo "-------------------------------"
echo "Removing redundant sequences..."

docker run --rm \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:x86_64-"${REFERENCE_PREPARATION_DOCKER_VERSION}" \
    bash -c "
        mmseqs createdb /data/tmp.filtered.fa /data/inputDB && \
        mmseqs clusthash /data/inputDB /data/resultDB --min-seq-id 1.0 && \
        mmseqs clust /data/inputDB /data/resultDB /data/clusterDB && \
        mmseqs createtsv /data/inputDB /data/inputDB /data/clusterDB /data/piRNA_cluster_result.tsv
    "

echo ""
echo "-------------------------------"
echo "Generating clustered FASTA file..."

docker run --rm \
    -v "${TMP_DIR}:/data" \
    -v "${CREATE_FASTA_MMSEQS2}:/scripts/create_fasta_mmseqs2.py" \
    ktrachtok/reference_preparation:x86_64-"${REFERENCE_PREPARATION_DOCKER_VERSION}" \
    python3 /scripts/create_fasta_mmseqs2.py \
        --fasta /data/tmp.filtered.fa \
        --mmseqs2_tsv /data/piRNA_cluster_result.tsv \
        --output /data/piRNA_db_custom.fa \
        --merge_headers

# Map piRNA to genome
# blat $GENOME $OUTPUT_DIR/piRNA_db_custom.fa -t=dna -q=rna -maxIntron=0 -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=100 -noTrimA -out=psl $OUTPUT_DIR/piRNA_db_custom_genomeMap.psl
echo ""
echo "-----------------------------------------"
echo "Alignment of piRNA sequences to genome..."
docker run --rm \
    -v "${TMP_DIR}:/data" \
    -v $GENOME:/data/genome.fa \
    ktrachtok/reference_preparation:x86_64-"${REFERENCE_PREPARATION_DOCKER_VERSION}" \
    ./blat /data/genome.fa /data/piRNA_db_custom.fa \
    -t=dna -q=rna -maxIntron=0 -stepSize=5 -repMatch=2253 \
    -minScore=20 -minIdentity=100 -noTrimA -out=psl \
    /data/piRNA_db_custom_genomeMap.psl

# Convert PSL to BED12
# Rscript $PSL2BED $OUTPUT_DIR/piRNA_db_custom_genomeMap.psl $OUTPUT_DIR/piRNA_db_custom_genomeMap.bed $OUTPUT_DIR
echo ""
echo "------------------------"
echo "Converting PSL to BED..."
docker run --rm \
    -v ${PSL2BED}:/scripts/psl2bed.r \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:x86_64-"${REFERENCE_PREPARATION_DOCKER_VERSION}" \
    Rscript /scripts/psl2bed.r /data/piRNA_db_custom_genomeMap.psl /data/piRNA_db_custom_genomeMap.bed

# Convert BED12 to GTF
# python $BED2GTF -i ${OUTPUT_DIR}/piRNA_db_custom_genomeMap.bed -o ${OUTPUT_DIR}/piRNA_db_custom_genomeMap.gtf --gene_feature --gene_biotype piRNA

# Convert BED12 to GTF
echo ""
echo "------------------------"
echo "Converting BED to GTF..."
docker run --rm \
    -v ${BED2GTF}:/scripts/bed2gtf.py \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:x86_64-"${REFERENCE_PREPARATION_DOCKER_VERSION}" \
    python3 /scripts/bed2gtf.py -i /data/piRNA_db_custom_genomeMap.bed -o /data/piRNA_db_custom_genomeMap.gtf --gene_feature --gene_biotype piRNA --source piRNA_custom_db

# STATISTICS OF DATABASES ####
echo ""
echo "------------------------"
echo "Generating statistics..."
# Get number of sequences in every source database
echo "$(basename $DB1),`grep -c ">" $DB1`" > $OUTPUT_DIR/pirna_databases.csv
echo "$(basename $DB2),`grep -c ">" $DB2`" >> $OUTPUT_DIR/pirna_databases.csv
echo "$(basename $DB3),`grep -c ">" $DB3`" >> $OUTPUT_DIR/pirna_databases.csv
echo "$(basename $DB4),`grep -c ">" $DB4`" >> $OUTPUT_DIR/pirna_databases.csv

# Get number of sequences after cleaning
echo "DB_Nclean,`grep -c ">" ${TMP_DIR}/tmp.filtered.fa`" >> $OUTPUT_DIR/pirna_databases.csv

# Get number of sequences after removing redundant piRNA sequences
echo "DB_uniqueSeq,`grep -c ">" $TMP_DIR/piRNA_db_custom.fa`" >> $OUTPUT_DIR/pirna_databases.csv

# Get number of piRNA that aligned to genome (with custom BLAT settings)
echo "DB_piRNAMapAll,$((`cut -f10 $TMP_DIR/piRNA_db_custom_genomeMap.psl | sort | uniq | wc -l`-5))" >> $OUTPUT_DIR/pirna_databases.csv

# Get number of piRNA that aligned to genome without any softclipping (full alignment)
echo "DB_piRNAMapNoClipping,`cut -f4 $TMP_DIR/piRNA_db_custom_genomeMap.bed | sort | uniq | wc -l`" >> $OUTPUT_DIR/pirna_databases.csv

# Length distribution of final piRNA database
echo "piRNA_length piRNA_number" > $OUTPUT_DIR/pirna_database_lenDist.csv
cat $TMP_DIR/piRNA_db_custom.fa | awk 'NR%2 == 0 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' >> $OUTPUT_DIR/pirna_database_lenDist.csv

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
' "$TMP_DIR/piRNA_db_custom.fa"
)

# Print the results
echo "Number of piRNA sequences: $NUM_SEQ"
echo "Minimum length: $MIN_LENGTH"
echo "Maximum length: $MAX_LENGTH"
echo "Average length: $AVG_LENGTH"

cp ${TMP_DIR}/piRNA_db_custom.fa ${OUTPUT_DIR}/piRNA_db_custom.fa
cp ${TMP_DIR}/piRNA_cluster_result.tsv ${OUTPUT_DIR}/piRNA_cluster_result.tsv
cp ${TMP_DIR}/piRNA_db_custom_genomeMap.gtf ${OUTPUT_DIR}/piRNA_db_custom_genomeMap.gtf
cp ${TMP_DIR}/piRNA_db_custom_genomeMap.bed ${OUTPUT_DIR}/piRNA_db_custom_genomeMap.bed
cp ${TMP_DIR}/piRNA_db_custom_genomeMap.psl ${OUTPUT_DIR}/piRNA_db_custom_genomeMap.psl

# Cleaning
rm -rf $TMP_DIR
