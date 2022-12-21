import pytest
from dna_sdr.sequence.dna_utilits import DNA
from dna_sdr.sequence.seq_utilits import sequence_comparison


t1 = DNA("TCT CA CAATC CA TCT")
base_loc_lst = [loc + 1 for loc in range(t1.base_count)]
zero_lst = [0] * t1.base_count
zero_dict = dict(zip(base_loc_lst, zero_lst))


@pytest.mark.parametrize("seq1, seq2, expected", [(t1, t1, zero_dict)])
def test_sequence_comparioson(seq1, seq2, expected):
    assert sequence_comparison(seq1, seq2) == expected
