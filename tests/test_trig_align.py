import pytest
from dna_sdr.sequence.seq_utilits import trig_aligment


# TODO: Update the test case
# TODO: Simplify the test file

test_lst = [
    ("CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA", (7, 5, 0, [], [])),
    ("CA TAACA CA TCT CA CATTC CA TCT CA CCACC CA", (7, 5, 1, [17], ["a-t"])),
    (
        "CA TAACA CA TCT CA CAATC CA TCT CA CCGGC CA",
        (7, 5, 2, [29, 30], ["a-g", "c-g"]),
    ),
    (
        "CA TAACA CA TCT CA CAATC CA TCT CA CAGGC CA",
        (7, 5, 3, [28, 29, 30], ["c-a", "a-g", "c-g"]),
    ),
    (
        "CA TAACA CA TCT CA CATTC CA TCT CA CCTCC CA",
        (7, 5, 2, [17, 29], ["a-t", "a-t"]),
    ),
    (
        "CA TAACA CA TCT CA CGTCC CA TCT CA CGTGC CA",
        (7, 5, 6, [16, 17, 18, 28, 29, 30], ["a-g", "a-t", "t-c", "c-g", "a-t", "c-g"]),
    ),
    ("CA CATCT CA CAATC CA TCT CA CCACC CA", (7, 0, 0, [], [])),
    ("CATCT CA CAATC CA TCT CA CCACC CA", (5, 0, 0, [], [])),
    ("CA TAACA TG TCT CA CAATC CA TCTCA CCACC CA", (5, 5, 2, [8, 9], ["c-t", "a-g"])),
    ("CA TAACA CT TCT CA CAGTC CA TCTCA CCACC CA", (6, 5, 2, [9, 17], ["a-t", "a-g"])),
    ("CA TAAAA CA TCT CA CAATC GA TCT CA CCACC CA", (6, 5, 2, [6, 20], ["c-a", "c-g"])),
    (
        "CA TAACT TG TCT CA CGTCC CA TCT CA CCACC CA",
        (4, 5, 6, [7, 8, 9, 16, 17, 18], ["a-t", "c-t", "a-g", "a-g", "a-t", "t-c"]),
    ),
    ("CAACT CA CAATC CA TCT CA CCACC CA", (5, 0, 1, [10], ["t-a"])),
    ("CAATT CA CAATC CA TCT CA CCACC CA", (5, 0, 2, [10, 11], ["t-a", "c-t"])),
    ("CA TAAAA CA TCT CT CAATC CA TCT CA CCACC CA", (6, 5, 2, [6, 14], ["c-a", "a-t"])),
    ("CA CATCT FA CAA TC CATCT CA CCACC CA", (ValueError)),
]


@pytest.mark.parametrize("trigger, expected", test_lst)
def test_strand_align(Comp_Trigger, trigger, expected):
    try:
        assert trig_aligment(Comp_Trigger, trigger) == expected
    except:
        with pytest.raises(TypeError):
            trig_aligment(Comp_Trigger, trigger)
