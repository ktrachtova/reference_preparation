# BED2GTF script
# Author Karolina Trachtova
# Takes BED12 formatted file and converts it to the GTF file with genes/transcripts/exons
# RUN:  python bed2gtf.py -i ~/smallRNA-Seq/pipeline_dev/databases/tRNA/genome_mapping/23_09_2022/tRNA_db_custom_gtRNAdb_genomeMap.bed -o ~/smallRNA-Seq/pipeline_dev/databases/tRNA/genome_mapping/23_09_2022/tRNA_db_custom_gtRNAdb_genomeMap.gtf --gene_feature --gene_biotype tRNA

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, action='store', dest='input', help="Input BED12 file")
    parser.add_argument("-o", "--output", action='store', dest='output', help="Output GTF file")
    parser.add_argument("-g", "--gene_feature", action='store_true', dest='gene_feature',
                        help="Should gene feature be reported as well? Will be identical to the transcript feature")
    parser.add_argument("-b", "--gene_biotype", type=str, action='store', dest='biotype', default="NA",
                        help="Gene biotype [default 'NA']")
    parser.add_argument("-r", "--add_group", action='store_true', dest='group')
    parser.add_argument("-s", "--source", action='store_true', dest='source')

    global args
    args = parser.parse_args()

    run()


# Function that will write a feature (gene/transcript/exon) into a file of choice
# All the necessary information about the feature are given as function parameters
# Required for all: chromosome, name of feature, start, end, strand, type of feature (gene/transcript/exon)
# Optional: group number (is user wants it, always present but wrote into file per user request),
# exon_number (0 for gene/transcript, not written for these features)
def writeFeature(chrom, source, name_gene, start, end, strand, feature, current_location, exon_number, fout):
    fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
        chrom,
        source,
        feature,
        str(start),
        str(end),
        ".",
        strand,
        ".",
        'gene_id "' + name_gene + '_loc' + current_location + '"; '
        + 'gene_name "' + name_gene + '_loc' + current_location + '"; '
        + 'gene_type "' + args.biotype + '";' +
        (
            ' transcript_id "' + name_gene + '_loc' + current_location + '"; ' + 'transcript_name "' + name_gene + '_loc' + current_location + '"; ' + 'transcript_type "' + args.biotype + '";' if feature == "transcript" or feature == "exon" else '') +
        (' exon_number "' + str(exon_number) + '";' if feature == "exon" else '')) + "\n")


def run():
    location_per_rna = {} # dictionary where number of genomic locations into which specific rna aligned to is stored

    fout = open(args.output, "w")

    with open(args.input, "r") as bed:
        for line in bed:
            lsplit = line.split("\t")
            chrom = lsplit[0]
            name_gene = lsplit[3]
            strand = lsplit[5]
            start = int(lsplit[1]) + 1
            end = int(lsplit[2])
            exon_sizes = lsplit[10].split(",")
            exon_starts = lsplit[11].split(",")

            if name_gene in location_per_rna.keys():
                new_location = location_per_rna[name_gene] + 1
                location_per_rna[name_gene] = new_location
            else:
                new_location = 1
                location_per_rna[name_gene] = 1

            if args.gene_feature:
                writeFeature(chrom, args.source, name_gene, start, end, strand, "gene", str(new_location), 0, fout)

            # write 'transcript' feature
            writeFeature(chrom, args.source, name_gene, start, end, strand, "transcript", str(new_location), 0, fout)

            # write 'exon' feature
            for i in range(0, len(exon_sizes) - 1):
                if i == 0:
                    exon_start = start
                    exon_end = start + int(exon_sizes[0]) - 1
                else:
                    exon_start = start + int(exon_starts[i])
                    exon_end = start + int(exon_starts[i]) + int(exon_sizes[i]) - 1

                writeFeature(chrom, args.source, name_gene, exon_start, exon_end, strand, "exon", str(new_location), i + 1, fout)


if __name__ == "__main__":
    main()
