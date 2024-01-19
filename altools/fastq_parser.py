from typing import Optional
from altools.dna_sequence import NamedDNASequence
from altools.nucleotide import Nucleotide, str_to_nucleotide

class FastQParserState:
    """Represents the internal state of the FASTQ Parser state machine"""

class Header(FastQParserState):
    def __init__(self, name: str):
        self.name = name
    
    def next(self, line: str):
        return DNA(self.name, str_to_nucleotide(line.strip())), None

class DNA(Header):
    def __init__(self, name: str, dna: tuple[Nucleotide]):
        super().__init__(name)
        self.dna = dna

    def next(self, _: str):
        return Separator(self.name, self.dna), None

class Separator(DNA):
    def __init__(self, name: str, dna: tuple[Nucleotide]):
        super().__init__(name, dna)

    def next(self, _: str):
        if self.name is not None and self.dna is not None:
            fastq_read = NamedDNASequence(self.name, self.dna)
        else:
            fastq_read = None
        return Quality(self.name, self.dna), fastq_read

class Quality(Separator):
    def __init__(self, name: str, dna: tuple[Nucleotide]):
        super().__init__(name, dna)

    def next(self, line: str):
        return Header(line.strip()), None

def read_fastq(fastq_path) -> tuple[NamedDNASequence]:
    """Reads in a FASTQ file"""

    state = Quality(None, None)
    
    with open(fastq_path) as fp:
        while True:
            line = fp.readline()
            if not line:
                return
            state, fastq_read = state.next(line)
            if fastq_read:
                yield fastq_read