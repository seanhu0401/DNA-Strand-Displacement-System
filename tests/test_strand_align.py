import pytest
from dna_sdr.sequence.seq_utilits import strand_alignment


@pytest.fixture
def default_tb_mismatch():
    return (True, 7)


# TODO: Simplify the test file

no_tb_mismatch_test_lst = [
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA", (7, 5, 0, [])),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCACC CA", (7, 5, 1, [17])),
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCGGC CA", (7, 5, 2, [29, 30])),
    ("CA TAACA CA TCT CA CAATC CA TCT CA CAGGC CA", (7, 5, 3, [28, 29, 30])),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCTCC CA", (7, 5, 2, [17, 29])),
    (
        "CA TAACA CA TCT CA CGTCC CA TCT CA CGTGC CA",
        (7, 5, 6, [16, 17, 18, 28, 29, 30]),
    ),
    ("CA CATCT CA CAATC CA TCT CA CCACC CA", (7, 0, 0, [])),
    ("CATCT CA CAATC CA TCT CA CCACC CA", (5, 0, 0, [])),
    ("CA CATCT FA CAATC CA TCT CA CCACC CA", (ValueError)),
]

tb_mismatch_test_lst = [
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA", (7, 5, 0, [])),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCACC CA", (7, 5, 1, [17])),
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCGGC CA", (7, 5, 2, [29, 30])),
    ("CA TAACA CA TCT CA CAATC CA TCT CA CAGGC CA", (7, 5, 3, [28, 29, 30])),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCTCC CA", (7, 5, 2, [17, 29])),
    (
        "CA TAACA CA TCT CA CGTCC CA TCT CA CGTGC CA",
        (7, 5, 6, [16, 17, 18, 28, 29, 30]),
    ),
    ("CA CATCT CA CAATC CA TCT CA CCACC CA", (7, 0, 0, [])),
    ("CA TAACA TG TCT CA CAATC CA TCTCA CCACC CA", (7, 5, 2, [8, 9])),
    ("CA TAACA CT TCT CA CAGTC CA TCTCA CCACC CA", (7, 5, 2, [9, 17])),
    ("CA CATCT FA CAA TC CATCT CA CCACC CA", (ValueError)),
]

tb_mismatch_test_not_defualt_lst = [
    ("CAACT CA CAA TC CATCT CA CCACC CA", (True, 5), (5, 0, 1, [3])),
    ("CAATT CA CAA TC CATCT CA CCACC CA", (True, 5), (5, 0, 2, [3, 4])),
    ("TC CAATT CA CAA TC CATCT CA CCACC CA", (True, 5), (5, 2, 2, [5, 6])),
    ("CA CATCT FA CAA TC CATCT CA CCACC CA", (True, 5), (ValueError)),
]


@pytest.mark.parametrize("trigger, expected", no_tb_mismatch_test_lst)
def test_strand_align_no_tb_mismatch(Incumb_Seq, Base_Seq, trigger, expected):
    try:
        assert strand_alignment(Base_Seq, Incumb_Seq, trigger) == expected
    except:
        with pytest.raises(TypeError):
            strand_alignment(Base_Seq, Incumb_Seq, trigger)


@pytest.mark.parametrize("trigger, expected", tb_mismatch_test_lst)
def test_strand_align_tb_mismatch(
    Incumb_Seq, Base_Seq, trigger, default_tb_mismatch, expected
):
    try:
        assert (
            strand_alignment(Base_Seq, Incumb_Seq, trigger, *default_tb_mismatch)
            == expected
        )
    except:
        with pytest.raises(TypeError):
            strand_alignment(Base_Seq, Incumb_Seq, trigger, *default_tb_mismatch)


@pytest.mark.parametrize(
    "trigger, tb_mismatch_option, expected", tb_mismatch_test_not_defualt_lst
)
def test_strand_align_tb_mismatch_not_default(
    Incumb_Seq, Base_Seq, trigger, tb_mismatch_option, expected
):
    try:
        assert (
            strand_alignment(Base_Seq, Incumb_Seq, trigger, *tb_mismatch_option)
            == expected
        )
    except:
        with pytest.raises(TypeError):
            strand_alignment(Base_Seq, Incumb_Seq, trigger, *tb_mismatch_option)
