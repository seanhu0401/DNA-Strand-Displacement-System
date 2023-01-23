import pytest
from dna_sdr.sequence.dna_utilits import DNA
from dna_sdr.sequence.seq_utilits import strand_alignment


@pytest.fixture
def tb_mismatch():
    return (True, 7)


no_tb_mismatch_test_lst = [
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA", (0, [], 5, 7)),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCACC CA", (1, [17], 5, 7)),
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCGGC CA", (2, [29, 30], 5, 7)),
    ("CA TAACA CA TCT CA CAATC CA TCT CA CAGGC CA", (3, [28, 29, 30], 5, 7)),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCTCC CA", (2, [17, 29], 5, 7)),
    (
        "CA TAACA CA TCT CA CGTCC CA TCT CA CGTGC CA",
        (6, [16, 17, 18, 28, 29, 30], 5, 7),
    ),
    ("CA CATCT CA CAATC CA TCT CA CCACC CA", (0, [], 0, 7)),
    ("CA CATCT FA CAATC CA TCT CA CCACC CA", (ValueError)),
]

tb_mismatch_test_lst = [
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA", (0, [], 5, 7)),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCACC CA", (1, [17], 5, 7)),
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCGGC CA", (2, [29, 30], 5, 7)),
    ("CA TAACA CA TCT CA CAATC CA TCT CA CAGGC CA", (3, [28, 29, 30], 5, 7)),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCTCC CA", (2, [17, 29], 5, 7)),
    (
        "CA TAACA CA TCT CA CGTCC CA TCT CA CGTGC CA",
        (6, [16, 17, 18, 28, 29, 30], 5, 7),
    ),
    ("CA CATCT CA CAA TC CATCT CA CCACC CA", (0, [], 0, 7)),
    ("CA TAACA TG TCT CA CAATC CA TCTCA CCACC CA", (2, [8, 9], 5, 7)),
    ("CA TAACA CT TCT CA CAGTC CA TCTCA CCACC CA", (2, [9, 17], 5, 7)),
    ("CA CATCT FA CAA TC CATCT CA CCACC CA", (ValueError)),
]


@pytest.mark.parametrize("trigger, expected", no_tb_mismatch_test_lst)
def test_strand_align_no_tb_mismatch(Incumb_Seq, Base_Seq, trigger, expected):
    try:
        assert strand_alignment(Base_Seq, Incumb_Seq, trigger) == expected
    except:
        with pytest.raises(TypeError):
            strand_alignment(Base_Seq, Incumb_Seq, trigger)


@pytest.mark.parametrize("trigger, expected", tb_mismatch_test_lst)
def test_strand_align_tb_mismatch(Incumb_Seq, Base_Seq, trigger, tb_mismatch, expected):
    try:
        assert strand_alignment(Base_Seq, Incumb_Seq, trigger, *tb_mismatch) == expected
    except:
        with pytest.raises(TypeError):
            strand_alignment(Base_Seq, Incumb_Seq, trigger, *tb_mismatch)
