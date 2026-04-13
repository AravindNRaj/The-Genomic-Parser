import sys
from collections import defaultdict
from itertools import product
import json

# IUPAC DNA ALPHABET
IUPAC_DNA = set("ATCGRYSWKMBDHVN")

# STANDARD BASE COMPLEMENTS
COMP = str.maketrans("ATCG", "TAGC")

# CODON TABLE (RNA → Amino Acid)
CODON_TABLE = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "AGU":"S","AGC":"S",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "UAU":"Y","UAC":"Y",
    "CAU":"H","CAC":"H",
    "CAA":"Q","CAG":"Q",
    "AAU":"N","AAC":"N",
    "AAA":"K","AAG":"K",
    "GAU":"D","GAC":"D",
    "GAA":"E","GAG":"E",
    "UGU":"C","UGC":"C",
    "UGG":"W",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGA":"R","AGG":"R",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
    "UAA":"STOP","UAG":"STOP","UGA":"STOP"
}

# VALIDATION
def validate(seq):
    invalid = set(seq) - IUPAC_DNA
    return len(invalid) == 0


# IUPAC MAP
IUPAC_MAP = {
    "A":["A"], "T":["T"], "C":["C"], "G":["G"],
    "R":["A","G"], "Y":["C","T"],
    "S":["G","C"], "W":["A","T"],
    "K":["G","T"], "M":["A","C"],
    "B":["C","G","T"], "D":["A","G","T"],
    "H":["A","C","T"], "V":["A","C","G"],
    "N":["A","T","C","G"]
}

def expand_codon(codon):
    bases = [IUPAC_MAP.get(b, [b]) for b in codon]
    return ["".join(p) for p in product(*bases)]

def translate_codon(codon):
    possibilities = expand_codon(codon)

    aas = set()
    for c in possibilities:
        aa = CODON_TABLE.get(c)
        if aa and aa != "STOP":
            aas.add(aa)

    if len(aas) == 1:
        return list(aas)[0]
    elif len(aas) > 1:
        return "X"
    else:
        return ""


# BIO FUNCTIONS
def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) * 100 if seq else 0

def transcribe(seq):
    return seq.replace("T", "U")

def reverse_complement(seq):
    return seq.translate(COMP)[::-1]


# ORF FINDER
def find_orfs(dna, strand="+"):

    seq = dna if strand == "+" else reverse_complement(dna)
    rna = transcribe(seq)

    orfs = []

    for frame in range(3):

        i = frame

        while i < len(rna) - 2:

            if rna[i:i+3] == "AUG":

                start = i
                protein = ""
                j = i

                while j < len(rna) - 2:

                    codon = rna[j:j+3]
                    aa = CODON_TABLE.get(codon, "")

                    if aa == "STOP":
                        end = j + 3

                        if len(protein) >= 5:

                            if strand == "+":
                                g_start = start
                                g_end = end
                                frame_label = f"+{frame+1}"
                            else:
                                g_start = len(dna) - end
                                g_end = len(dna) - start
                                frame_label = f"-{frame+1}"

                            score = len(protein) * (gc_content(dna[g_start:g_end]) + 1)

                            orfs.append({
                                "strand": strand,
                                "frame": frame_label,
                                "start": g_start,
                                "end": g_end,
                                "length_aa": len(protein),
                                "gc_score": round(gc_content(dna[g_start:g_end]), 2),
                                "score": round(score, 2),
                                "protein": protein
                            })

                        break

                    if aa == "":
                        break

                    protein += aa
                    j += 3

                i = j

            else:
                i += 3

    return orfs


# EXPORT GFF3
def export_gff3(orfs, filename):
    with open(filename, "w") as f:
        f.write("##gff-version 3\n")

        for i, o in enumerate(orfs, 1):
            attrs = (
                f"ID=orf{i};"
                f"length_aa={o['length_aa']};"
                f"score={o['score']}"
            )

            f.write(
                f"seq\tGenomicTool\tORF\t{o['start']}\t{o['end']}\t"
                f"{o['score']}\t{o['strand']}\t.\t{attrs}\n"
            )


# JSON EXPORT
def export_json(data, filename):
    with open(filename, "w") as f:
        json.dump(data, f, indent=4)


# FASTA PARSER
def read_fasta(filename):
    seqs = {}
    header = None
    seq = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq).upper()
                header = line[1:]
                seq = []
            else:
                seq.append(line)

        if header:
            seqs[header] = "".join(seq).upper()

    return seqs


# ANALYSIS
def analyze(header, seq):

    print("\n====================")
    print(header)
    print("====================")

    if not validate(seq):
        print("Invalid sequence skipped.")
        return None

    print(f"Length: {len(seq)} bp")
    print(f"GC Content: {gc_content(seq):.2f}%")

    forward = find_orfs(seq, "+")
    reverse = find_orfs(seq, "-")

    all_orfs = forward + reverse
    all_orfs.sort(key=lambda x: x["score"], reverse=True)

    print(f"Total ORFs: {len(all_orfs)}")

    if all_orfs:
        export_gff3(all_orfs, f"{header}_output.gff3")
        export_json(all_orfs, f"{header}_output.json")

        print("\nTop ORFs:")
        for o in all_orfs[:3]:
            print(f"{o['frame']} | {o['length_aa']} aa | score {o['score']}")

    return {header: all_orfs}


# MAIN
def main():

    if len(sys.argv) != 2:
        print("Usage: python Genome.py <file.fasta>")
        return

    data = read_fasta(sys.argv[1])

    if not data:
        print("No sequences found.")
        return

    results = {}

    for h, s in data.items():
        res = analyze(h, s)
        if res:
            results.update(res)

    export_json(results, "ALL_RESULTS.json")


if __name__ == "__main__":
    main()