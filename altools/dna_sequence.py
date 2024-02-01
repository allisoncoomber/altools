from typing import Iterable

from altools.nucleotide import str_to_nucleotide, Nucleotide

class DNASequence:
    def __init__(self, dna: Iterable[Nucleotide]):
        self.dna = tuple(dna)

    def __eq__(self, other: "DNASequence"):
        if isinstance(other, self.__class__):
            return self.dna == other.dna
        else:
            return False
    
    def __getitem__(self, key: int) -> Nucleotide: 
        return self.dna[key]
    
    def __len__(self) -> int:
        return len(self.dna)
    
    def __str__(self) -> str:
        return "".join([nuc.value for nuc in self.dna])

    @classmethod
    def from_string(cls, dna: str):
        return cls(str_to_nucleotide(dna))
    
class NamedDNASequence(DNASequence):
    def __init__(self, name: str, dna: tuple[Nucleotide]):
        super().__init__(dna)
        self.name = name

    def __eq__(self, other: "NamedDNASequence"):
        return self.name == other.name and self.dna == other.dna