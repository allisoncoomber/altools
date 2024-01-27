from altools.align.alignment import GlobalAlignment
from altools.align import align_nw
from altools.dna_sequence import DNASequence


def test_align_nw():
    alignments = align_nw(
        DNASequence.from_string("GCATGCG"),
        DNASequence.from_string("GATTACA"),
        match_reward=1,
        mismatch_penalty=-1,
        gap_penalty=-1,
    )
    expected_alignments = [
        GlobalAlignment(
            DNASequence.from_string("GCATG-CG"), DNASequence.from_string("G-ATTACA")
        ),
        GlobalAlignment(
            DNASequence.from_string("GCA-TGCG"), DNASequence.from_string("G-ATTACA")
        ),
        GlobalAlignment(
            DNASequence.from_string("GCAT-GCG"), DNASequence.from_string("G-ATTACA")
        ),
    ]

    for alignment in alignments:
        assert alignment in expected_alignments
    assert len(alignments) == len(expected_alignments)
