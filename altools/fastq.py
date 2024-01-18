from enum import Enum, auto
from typing import Optional

from altools.nucleotide import str_to_nucleotide, Nucleotide


class FASTQRead:
    def __init__(self, name: str, dna: tuple[Nucleotide]):
        self.name = name
        self.dna = dna
    def __eq__(self, other: "FASTQRead"):
        return self.name == other.name and self.dna == other.dna