from altools.align.alignment import Alignment
from altools.align.nw_semi_global import semi_global_align_nw
from altools.dna_sequence import DNASequence


def test_semi_global_nw_gaps_on_both_ends():
    alignments = semi_global_align_nw(
        DNASequence.from_string("GGGATAGGG"), DNASequence.from_string("ATA")
    )
    expected = [
        Alignment(
            DNASequence.from_string("GGGATAGGG"), DNASequence.from_string("---ATA---")
        )
    ]
    assert alignments == expected

def test_semi_global_nw_gaps_on_left():
    alignments = semi_global_align_nw(
        DNASequence.from_string("GGGATAGGG"), DNASequence.from_string("ATAGGG")
    )
    expected = [
        Alignment(
            DNASequence.from_string("GGGATAGGG"), DNASequence.from_string("---ATAGGG")
        )
    ]
    assert alignments == expected

def test_semi_global_nw_gaps_on_right():
    alignments = semi_global_align_nw(
        DNASequence.from_string("CCCATAGGG"), DNASequence.from_string("CCCATA")
    )
    expected = [
        Alignment(
            DNASequence.from_string("CCCATAGGG"), DNASequence.from_string("CCCATA---")
        )
    ]
    assert alignments == expected