#!/bin/bash

set -exo pipefail

# Set up project directory
PROJECT_DIR=$(dirname $(pwd))
echo "Project directory: $PROJECT_DIR"
mkdir -p $PROJECT_DIR # folder for downloaded reference files

GTF=""  # Gencode GTF file, example: gencode.v47.primary_assembly.annotation.gtf

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --GTF) GTF="$2"; shift ;;
        --help)
            echo "Usage: $0 --GTF <file>"
            exit 0
            ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Validate GTF
if [ -z "$GTF" ]; then
    echo "Error: --GTF (Gencode GTF file) is required."
    exit 1
fi

# Decompress GTF if gzipped
if [[ "$GTF" == *.gz ]]; then
    gzip -d -c "$GTF" > "${GTF%.gz}"
    GTF="${GTF%.gz}"
fi
echo "Using Gencode GTF file: $GTF"

# Create temporary folder for intermediate files
TMP_DIR="$(pwd)/tmp"
mkdir -p $TMP_DIR

# Script for filtering
FILTER_GTF=$(pwd)/utils/filter_gtf.py

# Output directory for results
OUTPUT_DIR=${PROJECT_DIR}/reference_files/other_rna/$(date +'%Y_%m_%d')
echo "${OUTPUT_DIR}"
mkdir -p $OUTPUT_DIR

# Extracting mRNA and lncRNA GTF records
#######################################################################################
echo ""
echo "----------------------------"
echo "Filtering GTF file to remove snoRNA/scaRNA/rRNA/Mt-tRNA features.."
echo ""

# Extract only protein_coding gene_type
OUTPUT_FILE=$(basename $GTF)
python3 ${FILTER_GTF} $GTF ${TMP_DIR}/${OUTPUT_FILE%.gtf}.filtered.gtf

echo ""
echo "Done"

cp ${TMP_DIR}/${OUTPUT_FILE%.gtf}.filtered.gtf ${OUTPUT_DIR}/${OUTPUT_FILE%.gtf}.filtered.gtf

# Clean
rm -rf $TMP_DIR
