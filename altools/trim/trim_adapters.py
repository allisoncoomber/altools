from altools.dna_sequence import DNASequence
from altools.align.nw_semi_global import semi_global_align_nw
from altools.align.alignment import Alignment


def trim_adapters(read: DNASequence, adapters: list[DNASequence]) -> DNASequence:
    """ Removes adapter(s) from the beginning or end of a DNA sequence

    For example
    Read: AAAAAGTC
    Adapter: GTC
    
    Would trim off the GTC from the read returning
    AAAAA

    Also:
    Read: AAAGAC
    Adapter: GTC
    Would not trim off the GTC from the read returning
    AAA
    """
    semi_global_align_nw(read, adapters[0])


def cut_alignment(alignment: Alignment) -> DNASequence:
