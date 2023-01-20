import pytest
from dna_sdr.sequence.dna_utilits import DNA
from dna_sdr.sequence.seq_utilits import strand_alignment

incumb = "CA CAATC CA TCT CA CCACC CA"
b_strand = "TG GGTGG TG AGA TG GATTG TG AGA TG TG AGA CAT ACA GCG CCG ACC GTA"

# TODO: parametrize the test file potentially


def test_no_mismatch_upper_case():
    trigger = "CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger, 7
    )
    assert len(mismatch_loc) == 0
    assert overhang_bs == 5


def test_no_mismatch_lower_case():
    trigger = "CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA"
    lower_trig = trigger.lower()
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger, 7
    )
    assert len(mismatch_loc) == 0
    assert overhang_bs == 5


def test_one_mismatch_no_gap():
    trigger = "CA TAACA CA TCT CA CATTC CA TCT CA CCACC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger, 7
    )
    assert len(mismatch_loc) == 1
    assert overhang_bs == 5


def test_two_mismatch_no_gap():
    trigger = "CA TAACA CA TCT CA CAATC CA TCT CA CCGGC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger, 7
    )
    assert len(mismatch_loc) == 2
    assert overhang_bs == 5


def test_three_mismatch_no_gap():
    trigger = "CA TAACA CA TCT CA CAATC CA TCT CA CAGGC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger, 7
    )
    assert len(mismatch_loc) == 3
    assert overhang_bs == 5


def test_two_mismatch_one_gap():
    trigger = "CA TAACA CA TCT CA CATTC CA TCT CA CCTCC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger, 7
    )
    assert len(mismatch_loc) == 2
    assert overhang_bs == 5


def test_six_mismatch_one_gap():
    trigger = "CA TAACA CA TCT CA CGTCC CA TCT CA CGTGC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger, 7
    )
    assert len(mismatch_loc) == 6
    assert overhang_bs == 5


def test_three_mismatch_three_gap():
    trigger = "CA TAACA CA TCT TA CAGTC CG TCT CA CCACC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger, 7
    )
    assert len(mismatch_loc) == 3
    assert overhang_bs == 5


def test_no_mismatch_no_gap_no_overhang():
    trigger = "CA CA TCT CA CAATC CA TCT CA CCACC CA "
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger, 7
    )
    assert len(mismatch_loc) == 0
    assert overhang_bs == 0


def test_not_fully_DNA():
    trigger = "CA CA TCT FA CAATC CA TCT CA CCACC CA"
    with pytest.raises(TypeError) as exc_info:
        mismatch_loc, overhang_bs, toehold_base = strand_alignment(
            incumb, b_strand, trigger, 7
        )
    exception_raised = exc_info.type
    assert exception_raised is TypeError


def test_no_mismatch_upper_case():
    trigger = "CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger
    )
    assert len(mismatch_loc) == 0
    assert overhang_bs == 5
    assert toehold_base == 7


def test_no_mismatch_lower_case():
    trigger = "CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA"
    lower_trig = trigger.lower()
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger
    )
    assert len(mismatch_loc) == 0
    assert overhang_bs == 5
    assert toehold_base == 7


def test_one_mismatch_no_gap():
    trigger = "CA TAACA CA TCT CA CATTC CA TCT CA CCACC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger
    )
    assert len(mismatch_loc) == 1
    assert overhang_bs == 5
    assert toehold_base == 7


def test_two_mismatch_no_gap():
    trigger = "CA TAACA CA TCT CA CAATC CA TCT CA CCGGC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger
    )
    assert len(mismatch_loc) == 2
    assert overhang_bs == 5
    assert toehold_base == 7


def test_three_mismatch_no_gap():
    trigger = "CA TAACA CA TCT CA CAATC CA TCT CA CAGGC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger
    )
    assert len(mismatch_loc) == 3
    assert overhang_bs == 5
    assert toehold_base == 7


def test_two_mismatch_one_gap():
    trigger = "CA TAACA CA TCT CA CATTC CA TCT CA CCTCC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger
    )
    assert len(mismatch_loc) == 2
    assert overhang_bs == 5
    assert toehold_base == 7


def test_six_mismatch_one_gap():
    trigger = "CA TAACA CA TCT CA CGTCC CA TCT CA CGTGC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger
    )
    assert len(mismatch_loc) == 6
    assert overhang_bs == 5
    assert toehold_base == 7


def test_three_mismatch_three_gap():
    trigger = "CA TAACA CA TCT TA CAGTC CG TCT CA CCACC CA"
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger
    )
    assert len(mismatch_loc) == 3
    assert overhang_bs == 5
    assert toehold_base == 7


def test_no_mismatch_no_gap_no_overhang():
    trigger = "CA CA TCT CA CAATC CA TCT CA CCACC CA "
    mismatch_loc, overhang_bs, toehold_base = strand_alignment(
        incumb, b_strand, trigger
    )
    assert len(mismatch_loc) == 0
    assert overhang_bs == 0
    assert toehold_base == 7


def test_not_fully_DNA():
    trigger = "CA CA TCT FA CAATC CA TCT CA CCACC CA"
    with pytest.raises(TypeError) as exc_info:
        mismatch_loc, overhang_bs, toehold_base = strand_alignment(
            incumb, b_strand, trigger
        )
    exception_raised = exc_info.type
    assert exception_raised is TypeError
