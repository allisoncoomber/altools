from altools.fastq_parser import read_fastq
from altools.dna_sequence import NamedDNASequence
from altools.nucleotide import Nucleotide


def test_read_fastq_one_read(tmp_path):
    (tmp_path / "sample.fastq").write_text("@read_name\n" "GATTACA\n" "+\n" "FFFFFFF")

    expected = (
        NamedDNASequence(
            "@read_name",
            (
                Nucleotide.G,
                Nucleotide.A,
                Nucleotide.T,
                Nucleotide.T,
                Nucleotide.A,
                Nucleotide.C,
                Nucleotide.A,
            ),
        ),
    )
    actual = tuple(read_fastq(tmp_path / "sample.fastq"))
    assert actual == expected


def test_read_fastq_two_reads(tmp_path):
    (tmp_path / "sample.fastq").write_text(
        "@read_name\n"
        "GATTACA\n"
        "+\n"
        "FFFFFFF\n"
        "@read2_name\n"
        "GATTACA\n"
        "+\n"
        "FFFFFFF"
    )

    expected = (
        NamedDNASequence(
            "@read_name",
            (
                Nucleotide.G,
                Nucleotide.A,
                Nucleotide.T,
                Nucleotide.T,
                Nucleotide.A,
                Nucleotide.C,
                Nucleotide.A,
            ),
        ),
        NamedDNASequence(
            "@read2_name",
            (
                Nucleotide.G,
                Nucleotide.A,
                Nucleotide.T,
                Nucleotide.T,
                Nucleotide.A,
                Nucleotide.C,
                Nucleotide.A,
            ),
        ),
    )
    actual = tuple(read_fastq(tmp_path / "sample.fastq"))
    assert actual == expected
