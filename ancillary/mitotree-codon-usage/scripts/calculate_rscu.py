#!/usr/bin/env python3
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict
import pandas as pd

# NCBI translation tables relevant to mitochondrial work.
# Stop codons are marked '_' and excluded from RSCU calculation.
GENETIC_CODES = {
    1: {  # Standard (nuclear)
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'_','TAG':'_',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'_','TGG':'W',
        'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    },
    2: {  # Vertebrate Mitochondrial
        # Changes vs Table 1: ATA=Met, TGA=Trp, AGA=Stop, AGG=Stop
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'M','ATG':'M',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'_','TAG':'_',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'_','AGG':'_',
        'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    },
    5: {  # Invertebrate Mitochondrial
        # Changes vs Table 1: ATA=Met, TGA=Trp, AGA=Ser, AGG=Ser
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'M','ATG':'M',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'_','TAG':'_',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'S','AGG':'S',
        'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    },
}


def calculate_rscu(fasta_file, genetic_code_table):
    if genetic_code_table not in GENETIC_CODES:
        sys.stderr.write(
            f"Unsupported genetic code {genetic_code_table}. "
            f"Available: {sorted(GENETIC_CODES.keys())}\n"
        )
        sys.exit(1)

    genetic_code = GENETIC_CODES[genetic_code_table]
    codon_counts = defaultdict(int)

    # Initialise counts for all sense codons to 0
    for codon, aa in genetic_code.items():
        if aa != '_':
            codon_counts[codon] = 0

    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq).upper()
            # Trim to multiple of 3
            remainder = len(seq) % 3
            if remainder != 0:
                seq = seq[:-remainder]

            for i in range(0, len(seq), 3):
                codon = seq[i:i+3]
                if len(codon) == 3 and codon in genetic_code and genetic_code[codon] != '_':
                    codon_counts[codon] += 1
    except Exception as e:
        sys.stderr.write(f"Error reading {fasta_file}: {e}\n")
        sys.exit(1)

    # Calculate RSCU: RSCU_ij = X_ij / (1/n_i * sum_j X_ij)
    results = []
    aa_to_codons = defaultdict(list)
    for codon, aa in genetic_code.items():
        if aa != '_':
            aa_to_codons[aa].append(codon)

    for aa, codons in aa_to_codons.items():
        total_count = sum(codon_counts[c] for c in codons)
        num_synonymous = len(codons)

        for codon in codons:
            count = codon_counts[codon]
            rscu = (count * num_synonymous) / total_count if total_count > 0 else 0.0
            results.append({
                'Amino_Acid': aa,
                'Codon': codon,
                'Count': count,
                'RSCU': rscu,
            })

    df = pd.DataFrame(results).sort_values(by=['Amino_Acid', 'Codon'])
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate Relative Synonymous Codon Usage (RSCU) from a CDS FASTA."
    )
    parser.add_argument("fasta_file", help="Input CDS FASTA file")
    parser.add_argument(
        "--genetic_code", type=int, default=2,
        help=(
            "NCBI translation table number (default: 2 = Vertebrate Mitochondrial). "
            "Supported: 1 (Standard/nuclear), 2 (Vertebrate Mitochondrial), "
            "5 (Invertebrate Mitochondrial)."
        ),
    )
    args = parser.parse_args()

    df = calculate_rscu(args.fasta_file, args.genetic_code)
    df.to_csv(sys.stdout, sep='\t', index=False)
