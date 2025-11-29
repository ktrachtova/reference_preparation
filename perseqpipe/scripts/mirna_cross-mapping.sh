#!/bin/bash

gencode_gtf=$1  # general GENCODE GTF file (needed to obtain miRNA gene coordinates)
sncrna_gtf=$2  # custom sncRNA GTF file for PerSeqPIPE

mkdir -p tmp/

echo "Input Gencode GTF: $gencode_gtf"
echo "Input PerSeqPIPE sncRNA GTF: $sncrna_gtf"
sncrna_gtf_version=$(echo "${sncrna_gtf##*_}" | sed 's/.gtf//')
echo "Input PerSeqPIPE sncRNA GTF version: $sncrna_gtf_version"

gencode_gtf_name=$(basename $gencode_gtf)
sncrna_gtf_name=$(basename $sncrna_gtf)

echo "Downloading docker image for bedtools..."
DOCKER_IMAGE="quay.io/biocontainers/bedtools"
DOCKER_TAG="2.27.1--h077b44d_9"
docker pull ${DOCKER_IMAGE}:${DOCKER_TAG}
DOCKER_IMAGE_ID=$(docker images -q $DOCKER_IMAGE)
echo "Docker image ID: ${DOCKER_IMAGE_ID}"
echo "Docker image downloaded"

echo ""
echo "Parsing Gencode GTF..."
echo "Extracting miRNA gene features."
grep -w 'gene_type "miRNA"' $gencode_gtf | grep -w 'gene' > ./tmp/${gencode_gtf_name%.gtf}.mirna_gene.gtf
echo "Example of parsed miRNA-gene file:"
head ./tmp/${gencode_gtf_name%.gtf}.mirna_gene.gtf
echo "Parsing done."

echo ""
echo "Parsing PerSeqPIPE sncRNA GTF..."
echo "Extracting only gene features."
grep -w 'gene' $sncrna_gtf > ./tmp/${sncrna_gtf_name%.gtf}_gene.gtf
echo "Example of parsed sncRNA-gene file:"
head ./tmp/${sncrna_gtf_name%.gtf}_gene.gtf
echo "Parsing done."

echo ""
echo "Running bedtools to identify overlaps:"
chmod -R 777 *
docker run --rm -v "$(pwd)":"$(pwd)" -w "$(pwd)" $DOCKER_IMAGE_ID bedtools intersect -a ./tmp/${gencode_gtf_name%.gtf}.mirna_gene.gtf -b ./tmp/${sncrna_gtf_name%.gtf}_gene.gtf -wo > ./tmp/mirna_sncrna_overlap.tsv
echo "Done."

echo ""
echo "Parsing bedtools results..."

awk -F'\t' '
BEGIN {
    OFS = "\t"
    # optional header:
    print "A_chr","A_start","A_end","A_gene_id","A_gene_name","A_gene_type", \
          "B_chr","B_start","B_end","B_gene_id","B_gene_name","B_gene_type", \
          "overlap"
}

# Extract value for key (gene_id/gene_name/gene_type) from GTF attribute field
function get_attr(field, key,    n, i, a, kv, v) {
    n = split(field, a, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^ +/, "", a[i])             # trim leading spaces
        if (index(a[i], key) == 1) {      # line starts with key
            split(a[i], kv, " ")
            v = kv[2]                     # value is in quotes
            gsub(/"/, "", v)              # remove quotes
            return v
        }
    }
    return "."
}

{
    # A attributes (col 9)
    a_gene_id   = get_attr($9,  "gene_id")
    a_gene_name = get_attr($9,  "gene_name")
    a_gene_type = get_attr($9,  "gene_type")

    # B attributes (col 18)
    b_gene_id   = get_attr($18, "gene_id")
    b_gene_name = get_attr($18, "gene_name")
    b_gene_type = get_attr($18, "gene_type")

    # print:
    #   A_chr A_start A_end A_gene_id A_gene_name A_gene_type
    #   B_chr B_start B_end B_gene_id B_gene_name B_gene_type
    #   overlap
    print $1, $4, $5, a_gene_id, a_gene_name, a_gene_type, \
          $10, $13, $14, b_gene_id, b_gene_name, b_gene_type, \
          $19
}
' ./tmp/mirna_sncrna_overlap.tsv > ./tmp/mirna_sncrna_overlap.clean.tsv_tmp

# keep header + sort body
(head -n 1 ./tmp/mirna_sncrna_overlap.clean.tsv_tmp && \
 tail -n +2 ./tmp/mirna_sncrna_overlap.clean.tsv_tmp | sort -t$'\t' -k12,12) \
 > ./tmp/mirna_sncrna_overlap_${sncrna_gtf_version}.tsv

# mv file to output folder
mkdir -p ../annotation_gtf/
mv ./tmp/mirna_sncrna_overlap_${sncrna_gtf_version}.tsv ../annotation_gtf/

rm -rf ./tmp/

echo "Parsing bedtools results done. Output: ./tmp/mirna_sncrna_overlap_${sncrna_gtf_version}.tsv"
