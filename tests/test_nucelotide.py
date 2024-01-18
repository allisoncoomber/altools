from altools.nucleotide import Nucleotide, str_to_nucleotide
import pytest

def test_str_to_nucleotide():
    assert str_to_nucleotide("ATCG") == (Nucleotide.A, Nucleotide.T, Nucleotide.C, Nucleotide.G)

def test_invalid_str_to_nucleotide():
    with pytest.raises(ValueError):
        str_to_nucleotide("EVIL")