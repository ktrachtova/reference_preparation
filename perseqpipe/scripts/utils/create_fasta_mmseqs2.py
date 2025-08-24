import argparse
from collections import defaultdict
from Bio import SeqIO
import sys

def extract_id(header: str) -> str:
    """Extract canonical identifier from FASTA header."""

    # URS identifiers
    if header.startswith("URS"):
        return header.split()[0]

    # Ensembl transcript IDs
    elif header.startswith("ENST"):
        return header.split("|")[0]

    # NCBI gi|...|ref|... or gi|...|gb|...
    elif header.startswith("gi|"):
        parts = header.split("|")
        if len(parts) >= 4:
            return parts[3]  # accession (e.g., XR_..., NR_..., M23820.1)
        else:
            sys.exit(f"ERROR: Invalid NCBI header: {header}\n"
                     f"Expected at least 4 fields separated by '|'")

    # tRNA entries
    elif "tRNA" in header:
        for part in header.split():
            if "tRNA" in part:
                return part.replace("Homo_sapiens_", "")
    
    # piRNA entries
    elif "piR-hsa" or "hsa-piR" in header:
        return header.split()[0]
    
    else:
        header_list = header.split()[0]
        if header_list == 0 or header_list == 1:
            print(f"WARNING: Header {header} is suspicious, check it is processed correctly!")
        else:
            return header.split()[0]

def main(fasta_file, tsv_file, output_file, merge_headers=False):
    # Load clusters from TSV
    clusters = defaultdict(list)
    with open(tsv_file) as f:
        for line in f:
            rep, member = line.strip().split()
            clusters[rep].append(member)

    # Load FASTA (keep original headers)
    seqs = list(SeqIO.parse(fasta_file, "fasta"))

    with open(output_file, "w") as out:
        for rep, members in clusters.items():
            # find representative record in FASTA
            rep_record = next((r for r in seqs if rep in r.description), None)
            if not rep_record:
                sys.exit(f"ERROR: Representative {rep} not found in FASTA")

            if merge_headers:
                cleaned_ids = []
                for m in members:
                    record = next((r for r in seqs if m in r.description), None)
                    if record:
                        cleaned_ids.append(extract_id(record.description))
                    else:
                        sys.exit(f"ERROR: Member {m} not found in FASTA! Are you using correct FASTA with TSV file from MMseqs2?")

                if not cleaned_ids:
                    continue
                header_id = "|".join(sorted(set(cleaned_ids)))
            else:
                header_id = extract_id(rep_record.description)

            out.write(f">{header_id}\n{str(rep_record.seq)}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge MMseqs2 clusters with custom headers")
    parser.add_argument("--fasta", required=True, help="Input FASTA file (original sequences)")
    parser.add_argument("--mmseqs2_tsv", required=True, help="TSV file from mmseqs createtsv")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--merge_headers", action="store_true",
                        help="Merge all member IDs into header (default: use only representative ID)")
    args = parser.parse_args()

    main(args.fasta, args.mmseqs2_tsv, args.output, merge_headers=args.merge_headers)
