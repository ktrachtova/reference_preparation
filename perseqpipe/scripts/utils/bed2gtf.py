#!/usr/bin/env python3

# BED2GTF script
# Author: Karolina Trachtova

import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Input BED12 file")
    parser.add_argument("-o", "--output", required=True, help="Output GTF file")
    parser.add_argument("-g", "--gene_feature", action='store_true',
                        help="Should gene feature be reported as well? Will be identical to the transcript feature")
    parser.add_argument("-b", "--gene_biotype", type=str, default="NA",
                        help="Gene biotype [default 'NA']")
    parser.add_argument("-s", "--source", type=str, default="BED2GTF",
                        help="Source field for GTF [default 'BED2GTF']")

    args = parser.parse_args()
    run(args)


def writeFeature(chrom, source, gene_id, gene_name, start, end, strand, feature,
                 transcript_id, exon_number, biotype, fout):
    """Write a single GTF feature line."""
    attributes = [
        f'gene_id "{gene_id}";',
        f'gene_name "{gene_name}";',
        f'gene_type "{biotype}";'
    ]
    if feature in ("transcript", "exon"):
        attributes += [
            f'transcript_id "{transcript_id}";',
            f'transcript_name "{transcript_id}";',
            f'transcript_type "{biotype}";'
        ]
    if feature == "exon":
        attributes.append(f'exon_number "{exon_number}";')

    fout.write("\t".join([
        chrom,
        source,
        feature,
        str(start),
        str(end),
        ".",
        strand,
        ".",
        " ".join(attributes)
    ]) + "\n")


def run(args):
    # First pass: count number of alignments per RNA
    counts = defaultdict(int)
    with open(args.input) as bed:
        for line in bed:
            if not line.strip():
                continue
            fields = line.split("\t")
            name_gene = fields[3]
            counts[name_gene] += 1

    # Second pass: write features
    location_per_rna = defaultdict(int)

    with open(args.input) as bed, open(args.output, "w") as fout:
        for line in bed:
            if not line.strip():
                continue
            lsplit = line.split("\t")
            chrom = lsplit[0]
            name_gene = lsplit[3]
            strand = lsplit[5]
            start = int(lsplit[1]) + 1
            end = int(lsplit[2])
            exon_sizes = lsplit[10].split(",")
            exon_starts = lsplit[11].split(",")

            # increment mapping counter
            location_per_rna[name_gene] += 1
            loc_index = location_per_rna[name_gene]

            # decide ID format
            if counts[name_gene] > 1:
                gene_id = f"{name_gene}_loc{loc_index}"
                transcript_id = gene_id
                gene_name = gene_id
            else:
                gene_id = name_gene
                transcript_id = name_gene
                gene_name = name_gene

            if args.gene_feature:
                writeFeature(chrom, args.source, gene_id, gene_name,
                             start, end, strand, "gene", transcript_id, 0,
                             args.gene_biotype, fout)

            # transcript
            writeFeature(chrom, args.source, gene_id, gene_name,
                         start, end, strand, "transcript", transcript_id, 0,
                         args.gene_biotype, fout)

            # exons
            for i in range(len(exon_sizes) - 1):  # last entry is often empty
                if i == 0:
                    exon_start = start
                    exon_end = start + int(exon_sizes[0]) - 1
                else:
                    exon_start = start + int(exon_starts[i])
                    exon_end = exon_start + int(exon_sizes[i]) - 1

                writeFeature(chrom, args.source, gene_id, gene_name,
                             exon_start, exon_end, strand, "exon", transcript_id,
                             i + 1, args.gene_biotype, fout)


if __name__ == "__main__":
    main()
