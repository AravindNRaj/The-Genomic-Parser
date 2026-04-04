import sys
def read_fasta(filename):
    """Read a FASTA file and return the DNA sequence as a string."""
    sequence = ""
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue  # skip header
            sequence += line.strip().upper()
    return sequence

def gc_content(seq):
    """Calculate GC-content as a percentage."""
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100

def transcribe_dna(seq):
    """Convert DNA sequence to RNA sequence."""
    return seq.replace("T", "U")

def translate_rna(rna_seq):
    """Translate RNA into a protein string using the standard genetic code."""
    codon_table = {
        # Start codon
        "AUG": "M",
        # Stop codons
        "UAA": "*", "UAG": "*", "UGA": "*",
        # Example codons (you'll need to complete all 64)
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUU": "I", "AUC": "I", "AUA": "I",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        # ... fill in the rest
    }
    protein = ""
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        protein += codon_table.get(codon, "?")  # unknown codons as '?'
    return protein

def main():
    if len(sys.argv) != 2:
        print("Usage: python dna_tool.py <file.fasta>")
        sys.exit(1)

    filename = sys.argv[1]
    dna_seq = read_fasta(filename)
    print(f"GC-content: {gc_content(dna_seq):.2f}%")
    
    rna_seq = transcribe_dna(dna_seq)
    print(f"RNA sequence: {rna_seq}")
    
    protein_seq = translate_rna(rna_seq)
    print(f"Protein sequence: {protein_seq}")

if __name__ == "__main__":
    main()