#!/usr/bin/env python3
import sys
from collections import defaultdict
from datetime import datetime

# Set of unwanted gene types to filter out
FILTER_TYPES = {"snoRNA", "scaRNA", "rRNA", "Mt_tRNA", "miRNA"}

def parse_attributes(attr_str):
    """Parse the attributes column of GTF into a dict."""
    attrs = {}
    for kv in attr_str.strip().split(';'):
        kv = kv.strip()
        if not kv:
            continue
        if ' ' not in kv:
            continue
        key, value = kv.split(' ', 1)
        attrs[key] = value.strip('"')
    return attrs

def read_gtf(input_path):
    """Read GTF file line by line."""
    with open(input_path, "r") as f:
        for line in f:
            yield line.rstrip("\n")

def main(input_gtf, output_gtf):
    lines = []
    gene_to_type = {}
    transcript_to_gene = {}
    gene_to_transcripts = defaultdict(set)
    gene_to_features = defaultdict(set)

    # First pass: collect all metadata
    for line in read_gtf(input_gtf):
        if line.startswith("#"):
            continue

        cols = line.split("\t")
        if len(cols) < 9:
            continue

        feature_type = cols[2]
        attributes = parse_attributes(cols[8])
        gene_id = attributes.get("gene_id")
        transcript_id = attributes.get("transcript_id")
        gene_type = attributes.get("gene_type")

        if gene_id:
            gene_to_type[gene_id] = gene_type
            gene_to_features[gene_id].add(feature_type)
        if transcript_id and gene_id:
            transcript_to_gene[transcript_id] = gene_id
            gene_to_transcripts[gene_id].add(transcript_id)

        lines.append((line, gene_id))

    # Determine which genes to keep
    kept_genes = {g for g, t in gene_to_type.items() if t not in FILTER_TYPES}

    filtered_lines = []
    for line, gene_id in lines:
        # Keep headers always
        if line.startswith("#"):
            filtered_lines.append(line)
            continue

        # Drop records for removed genes
        if gene_id not in kept_genes:
            continue

        filtered_lines.append(line)

    # Final consistency check: no orphan features allowed
    for line in filtered_lines:
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        attributes = parse_attributes(cols[8])
        gene_id = attributes.get("gene_id")

        # If a feature refers to a removed gene, that's an error
        if gene_id not in kept_genes:
            print("\n" + "="*80, file=sys.stderr)
            print(f"ERROR: Orphan record found!", file=sys.stderr)
            print(f"  Gene ID: {gene_id}", file=sys.stderr)
            print(f"  Full attributes: {cols[8]}", file=sys.stderr)
            print(f"  Offending line: {line}", file=sys.stderr)
            print("="*80 + "\n", file=sys.stderr)
            sys.exit(1)

    # Write the final filtered GTF
    with open(output_gtf, "w") as out:
        out.write("\n".join(filtered_lines) + "\n")

    print(f"Filtering completed successfully. Output saved to: {output_gtf}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.gtf> <output.gtf>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
