#!/bin/bash

# Set up project directory
PROJECT_DIR=$(dirname $(pwd))
echo "Project directory: $PROJECT_DIR"
mkdir -p $PROJECT_DIR # folder for downloaded reference files

DB1=""  # RNACentral FASTA
DB2=""  # Gencode FASTA
GENOME=""  # GRCh38.primary_assembly.genome.fa 

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --DB1) DB1="$2"; shift ;;
        --DB2) DB2="$2"; shift ;;
        --GENOME) GENOME="$2"; shift ;;
        --help)
            echo "Usage: $0 --DB1 <file> --DB2 <file> [--GENOME <file>]"
            exit 0
            ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Validate DB1 and DB2
if [ -z "$DB1" ]; then
    echo "Error: --DB1 (RNACentral FASTA) is required."
    exit 1
fi

if [ -z "$DB2" ]; then
    echo "Error: --DB2 (GENCODE FASTA) is required."
    exit 1
fi

# Decompress DB1 if gzipped
if [[ "$DB1" == *.gz ]]; then
    gzip -d -c "$DB1" > "${DB1%.gz}"
    DB1="${DB1%.gz}"
fi
echo "Using RNACentral FASTA file: $DB1"

# Decompress DB2 if gzipped
if [[ "$DB2" == *.gz ]]; then
    gzip -d -c "$DB2" > "${DB2%.gz}"
    DB2="${DB2%.gz}"
fi
echo "Using GENCODE FASTA file: $DB2"

# Handle GENOME download if not provided
if [ -z "$GENOME" ]; then
    echo "Genome FASTA file not provided, downloading..."
    mkdir -p "$PROJECT_DIR/databases"
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
    gzip -d GRCh38.primary_assembly.genome.fa.gz
    mv GRCh38.primary_assembly.genome.fa "$PROJECT_DIR/databases"
    GENOME="$PROJECT_DIR/databases/GRCh38.primary_assembly.genome.fa"
else
    if [[ "$GENOME" == *.gz ]]; then
        gzip -d -c "$GENOME" > "${GENOME%.gz}"
        GENOME="${GENOME%.gz}"
    fi
    echo "Using genome FASTA: $GENOME"
fi

# PSL to BED12 R conversion scripts
PSL2BED=$(pwd)/utils/psl2bed.r
BED2GTF=$(pwd)/utils/bed2gtf.py
CREATE_FASTA_MMSEQS2=$(pwd)/utils/create_fasta_mmseqs2.py

# Create temporary folder for intermediate files
TMP_DIR="$(pwd)/tmp"
mkdir -p $TMP_DIR

# Output directory for results
OUTPUT_DIR=${PROJECT_DIR}/reference_files/snoRNA/$(date +'%Y_%m_%d')
echo "${OUTPUT_DIR}"
mkdir -p $OUTPUT_DIR

# File with sequences to remove
# REMOVE_SEQ_FILE=${PROJECT_DIR}/databases/snoRNA/rnacentral_snorna_remove.txt

# Create snoRNA custom database ####
#######################################################################################
echo ""
echo "----------------------------"
echo "Cleaning snoRNA sequences..."

# FASTA multiline to one line
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' $DB1 > $TMP_DIR/db1_
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' $DB2 > $TMP_DIR/db2_

# Replace Us with Ts
awk '{ if ($0 !~ /^>/) gsub(/U/, "T"); print }' $TMP_DIR/db1_ > $TMP_DIR/db1_Ts_
# awk '{ if ($0 !~ /^>/) gsub(/U/, "T"); print }' $TMP_DIR/db2_ > $TMP_DIR/db2_Ts_ -> Gencode has Ts, so not needed

# Extract only snoRNA transcripts from the GENCODE FASTA file
grep -E -A 1 "snoRNA|scaRNA" $TMP_DIR/db2_ | grep -v "^--" > $TMP_DIR/db2_snoRNA_

# Merge both databases
cat $TMP_DIR/db1_Ts_ $TMP_DIR/db2_snoRNA_ > $TMP_DIR/tmp.fa

echo ""
echo "-------------------------------"
echo "Removing redundant sequences..."

docker run --rm \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:latest \
    bash -c "
        /MMseqs2/build/bin/mmseqs createdb /data/tmp.fa /data/inputDB && \
        /MMseqs2/build/bin/mmseqs clusthash /data/inputDB /data/resultDB --min-seq-id 1.0 && \
        /MMseqs2/build/bin/mmseqs clust /data/inputDB /data/resultDB /data/clusterDB && \
        /MMseqs2/build/bin/mmseqs createtsv /data/inputDB /data/inputDB /data/clusterDB /data/snoRNA_cluster_result.tsv
    "

echo ""
echo "-------------------------------"
echo "Generating clustered FASTA file..."

docker run --rm \
    -v "${TMP_DIR}:/data" \
    -v "${CREATE_FASTA_MMSEQS2}:/scripts/create_fasta_mmseqs2.py" \
    ktrachtok/reference_preparation:latest \
    python3 /scripts/create_fasta_mmseqs2.py \
        --fasta /data/tmp.fa \
        --mmseqs2_tsv /data/snoRNA_cluster_result.tsv \
        --output /data/snoRNA_db_custom.fa \
        --merge_headers

echo ""
echo "--------------------------------------"
echo "Aligning snoRNA sequences to genome..."

docker run --rm \
    -v "${TMP_DIR}:/data" \
    -v $GENOME:/data/genome.fa \
    ktrachtok/reference_preparation:latest \
    ./blat /data/genome.fa /data/snoRNA_db_custom.fa \
    -t=dna -q=rna -repMatch=1000 \
    -minScore=30 -minIdentity=100 -noTrimA -fine -out=psl \
    /data/snoRNA_db_custom_genomeMap.psl


# Convert PSL to BED12
# Rscript $PSL2BED ${OUTPUT_DIR}/snoRNA_db_custom_genomeMap.psl ${OUTPUT_DIR}/snoRNA_db_custom_genomeMap.bed $OUTPUT_DIR

echo ""
echo "---------------------"
echo "Converting PSL to BED"

docker run --rm \
    -v ${PSL2BED}:/scripts/psl2bed.r \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:latest \
    Rscript /scripts/psl2bed.r /data/snoRNA_db_custom_genomeMap.psl /data/snoRNA_db_custom_genomeMap.bed

# Convert BED12 to GTF
# python $BED2GTF -i ${OUTPUT_DIR}/snoRNA_db_custom_genomeMap.bed -o ${OUTPUT_DIR}/snoRNA_db_custom_genomeMap.gtf --gene_feature --gene_biotype snoRNA

echo ""
echo "------------------------"
echo "Converting BED to GTF..."

docker run --rm \
    -v ${BED2GTF}:/scripts/bed2gtf.py \
    -v "${TMP_DIR}:/data" \
    ktrachtok/reference_preparation:latest \
    python3 /scripts/bed2gtf.py -i /data/snoRNA_db_custom_genomeMap.bed -o /data/snoRNA_db_custom_genomeMap.gtf --gene_feature --gene_biotype snoRNA --source snoRNA_custom_db

# Calculate and print some statistics about created snoRNA database
echo ""
echo "-------------------------------"
echo "Calculating final statistics..."

# Calculate minimum, maximum, and average lengths
MIN_LENGTH=$(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' "${TMP_DIR}/snoRNA_db_custom.fa" | sort -n | head -n 1)
MAX_LENGTH=$(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' "${TMP_DIR}/snoRNA_db_custom.fa" | sort -n | tail -n 1)
AVG_LENGTH=$(awk '/^>/ {if (seqlen) {sum += seqlen; count++} seqlen=0; next} {seqlen += length($0)} END {print sum/count}' "${TMP_DIR}/snoRNA_db_custom.fa")
NUM_SEQ=`grep -c ">" "${TMP_DIR}/snoRNA_db_custom.fa"`

# Print the results
echo "Number of unique snoRNA: $NUM_SEQ"
echo "Minimum length: $MIN_LENGTH"
echo "Maximum length: $MAX_LENGTH"
echo "Average length: $AVG_LENGTH"

# Calculate number of snoRNA sequences in source databases as well as the final database
echo "$(basename $DB1),`grep -c ">" $DB1`" > ${OUTPUT_DIR}/snorna_databases.csv
echo "$(basename $DB2),`grep -c ">" ${TMP_DIR}/db2_snoRNA_`" >> ${OUTPUT_DIR}/snorna_databases.csv

echo "DB_allSeq, `grep -c ">" $TMP_DIR/tmp.fa`" >> ${OUTPUT_DIR}/snorna_databases.csv
echo "DB_uniqueSeq,`grep -c ">" ${TMP_DIR}/snoRNA_db_custom.fa`" >> ${OUTPUT_DIR}/snorna_databases.csv

# length distribution of final tRNA database
echo "snoRNA_length snoRNA_number" > ${OUTPUT_DIR}/snorna_database_lenDist.csv
cat ${TMP_DIR}/snoRNA_db_custom.fa | awk 'NR%2 == 0 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' >> ${OUTPUT_DIR}/snorna_database_lenDist.csv

echo "Done."

cp ${TMP_DIR}/snoRNA_db_custom.fa ${OUTPUT_DIR}/snoRNA_db_custom.fa
cp ${TMP_DIR}/snoRNA_cluster_result.tsv ${OUTPUT_DIR}/snoRNA_cluster_result.tsv
cp ${TMP_DIR}/snoRNA_db_custom_genomeMap.gtf ${OUTPUT_DIR}/snoRNA_db_custom_genomeMap.gtf
cp ${TMP_DIR}/snoRNA_db_custom_genomeMap.bed ${OUTPUT_DIR}/snoRNA_db_custom_genomeMap.bed
cp ${TMP_DIR}/snoRNA_db_custom_genomeMap.psl ${OUTPUT_DIR}/snoRNA_db_custom_genomeMap.psl

# Clean
rm -rf $TMP_DIR
