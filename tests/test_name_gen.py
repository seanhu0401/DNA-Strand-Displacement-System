import pytest
from dna_sdr.sequence.dna_utilits import DNA
from dna_sdr.sequence.seq_utilits import trig_alignment
from dna_sdr.sequence.seq_utilits import name_generation

test_lst = [
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA", "7tb_5'5ob"),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCACC CA", "7tb_5'5ob_17_1mb_a-t"),
    (
        "CA TAACA CA TCT CA CAATC CA TCT CA CAGGC CA",
        "7tb_5'5ob_28_3mb_('c-a', 'a-g', 'c-g')",
    ),
    (
        "CA TAACA CA TCT CA CGTCC CA TCT CA CGTGC CA",
        "7tb_5'5ob_16_3mb_('a-g', 'a-t', 't-c')_28_3mb_('c-g', 'a-t', 'c-g')",
    ),
]


@pytest.mark.parametrize("trig, expected", test_lst)
def test_seq_comparision(Comp_Trigger, trig, expected):
    info = trig_alignment(Comp_Trigger, trig)
    assert name_generation(*info[:2], *info[3:]) == expected
