from enum import Enum

class Nucleotide(Enum):
    A = "A"
    C = "C"
    G = "G"
    T = "T"
    GAP = "-"

def str_to_nucleotide(dna: str) -> tuple[Nucleotide]:
    return tuple(Nucleotide(x) for x in dna)