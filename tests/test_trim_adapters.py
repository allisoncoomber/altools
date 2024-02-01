
from altools.dna_sequence import DNASequence
from altools.trim.trim_adapters import trim_adapters


def test_trim_adapters_one():
    actual = trim_adapters(DNASequence.from_string("GATTACA"), [DNASequence.from_string("ACA")])
    expected = DNASequence.from_string("GATT")
    assert actual == expected