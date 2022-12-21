import pytest
from dna_sdr.sequence.dna_utilits import DNA
from dna_sdr.sequence.seq_utilits import sequence_comparison


@pytest.fixture
def t1():
    return DNA("TCT CA CAATC CA TCT")


@pytest.fixture
def base_dict_gen(t1):
    lst = [0] * t1.base_count
    base_loc_lst = [loc + 1 for loc in range(t1.base_count)]
    base_dict = dict(zip(base_loc_lst, lst))
    return base_dict


test_lst = [
    (DNA("TCT CA CAATC CA TCT"), {}),
    (DNA("TCT CA TAATC CA TCT"), {6: 1}),
    (DNA("TCT CA GAATC CA TCT"), {6: 2}),
    (DNA("TCT CA AAATC CA TCT"), {6: 4}),
    (DNA("TCT CA CACTC CA TCT"), {8: 3}),
    (DNA("TCT CA CCCTC CA TCT"), {7: 3, 8: 3}),
    (DNA("TCT GA CAATA CA TCT"), {4: 2, 10: 4}),
]


@pytest.mark.parametrize("t2, updates", test_lst)
def test_seq_comparision(t1, base_dict_gen, t2, updates):
    expected = base_dict_gen
    expected.update(updates)
    assert sequence_comparison(t1, t2) == expected
