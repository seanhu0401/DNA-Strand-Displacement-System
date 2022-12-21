import pytest
from dna_sdr.sequence.dna_utilits import DNA
from dna_sdr.sequence.seq_utilits import sequence_comparison


@pytest.fixture
def t1():
    return DNA("TCT CA CAATC CA TCT")


def test_same_sequence(t1):
    zero_lst = [0] * t1.base_count
    base_loc_lst = [loc + 1 for loc in range(t1.base_count)]
    expected = dict(zip(base_loc_lst, zero_lst))
    assert sequence_comparison(t1, t1) == expected
