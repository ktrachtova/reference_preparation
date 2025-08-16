import argparse
from collections import defaultdict
from Bio import SeqIO

def extract_id(header):
    """Extracts identifier based on header pattern."""
    if header.startswith("URS"):
        return header.split()[0]
    elif header.startswith("ENST"):
        return header.split("|")[0]
    elif "tRNA" in header:
        for part in header.split():
            if "tRNA" in part:
                return part.replace("Homo_sapiens_", "")
    else:
        return header.split()[0]

def main(fasta_file, tsv_file, output_file):
    # Load sequences
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Load clusters
    clusters = defaultdict(list)
    with open(tsv_file) as f:
        for line in f:
            rep, member = line.strip().split()
            clusters[rep].append(member)

    # Write output
    with open(output_file, "w") as out:
        for rep, members in clusters.items():
            ids = [extract_id(seqs[m].description) for m in members]
            merged_id = "|".join(sorted(set(ids)))
            seq = str(seqs[rep].seq)  # one-line sequence
            out.write(f">{merged_id}\n{seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge MMseqs2 clusters with custom headers")
    parser.add_argument("--fasta", required=True, help="Input FASTA file (original sequences)")
    parser.add_argument("--mmseqs2_tsv", required=True, help="TSV file from mmseqs createtsv")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    args = parser.parse_args()
    main(args.fasta, args.mmseqs2_tsv, args.output)
