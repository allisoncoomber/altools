from enum import Enum, auto
from typing import Optional

from altools.nucleotide import str_to_nucleotide, Nucleotide

class DNASequence:
    def __init__(self, dna: tuple[Nucleotide]):
        self.dna = dna

    def __eq__(self, other: "DNASequence"):
        return self.dna == other.dna
    
    def __getitem__(self, key: int) -> Nucleotide: 
        return self.dna[key]
    
    def __len__(self) -> int:
        return len(self.dna)

    @classmethod
    def from_string(cls, dna: str):
        return cls(str_to_nucleotide(dna))
    
class NamedDNASequence(DNASequence):
    def __init__(self, name: str, dna: tuple[Nucleotide]):
        super().__init__(dna)
        self.name = name

    def __eq__(self, other: "NamedDNASequence"):
        return self.name == other.name and self.dna == other.dna