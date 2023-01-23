from dna_sdr.sequence.dna_utilits import DNA

comp_dict = {"t": 1, "g": 2, "c": 3, "a": 4}
base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}

# *Think of better name for the method
# This function is currently used to compare to sequence, mainly the triggers, to identify the location(s) that has/have mismatch(es)
# They should be equal length in bases if there is no modification to the toehold + dangles - not applicable with those modification currently
def sequence_comparison(seq_1: str, seq_2: str):
    comparison_lst = []
    base_loc_lst = []
    counter = 0

    try:
        strand_1 = DNA(seq_1)
        strand_2 = DNA(seq_2)
    except ValueError:
        raise TypeError("The input is not a DNA type")

    if strand_1.base_count != strand_2.base_count:
        raise Exception("The two sequences are not equal length")

    else:
        while counter < strand_1.base_count:
            comp_val = 0
            if strand_1.seq[counter] == strand_2.seq[counter]:
                comparison_lst.append(comp_val)
            elif strand_1.seq[counter] != strand_2.seq[counter]:
                comp_val = comp_dict[strand_2.seq[counter]]
                comparison_lst.append(comp_val)
            base_loc_lst.append(counter + 1)
            counter += 1
        seq_comp_dict = dict(zip(base_loc_lst, comparison_lst))

    return seq_comp_dict


def strand_alignment_no_tb_mismatch(base_strand: DNA, incumbent: DNA, trigger: DNA):
    bstrand = base_strand.seq
    b_base = base_strand.base_count

    b_incumb = incumbent.base_count

    reverse_trig = trigger.reverse_seq()
    b_trig = trigger.base_count

    counter = mismatch_count = 0
    last_match_loc = 0
    tb_ob = b_trig - b_incumb  # number of bases for toehold + overhang
    mismatch_loc = list()

    for base in reverse_trig:
        comp_base = base_pairing[base]
        base_at_loc = bstrand[counter]

        if comp_base == base_at_loc and counter <= b_base:
            pass

        elif comp_base != base_at_loc and b_trig - counter > tb_ob:
            base_loc_5 = b_trig - counter
            mismatch_loc.append(base_loc_5)
            mismatch_count += 1

        elif (
            comp_base != base_at_loc
            and b_trig - counter <= tb_ob
            and counter > last_match_loc
        ):
            if last_match_loc < b_trig - counter:
                last_match_loc = counter

        counter += 1

    if last_match_loc != 0:
        toehold = last_match_loc - b_incumb
        overhang = tb_ob - toehold

    else:
        toehold = tb_ob
        overhang = 0

    mismatch_loc.sort()

    return mismatch_count, mismatch_loc, overhang, toehold


def strand_alignment_tb_mismatch(base_strand: DNA, incumbent: DNA, trigger: DNA, tb: 7):
    bstrand = base_strand.seq
    b_base = base_strand.base_count

    b_incumb = incumbent.base_count

    reverse_trig = trigger.reverse_seq()
    b_trig = trigger.base_count

    counter = mismatch = 0
    base_count = b_incumb + tb
    mismatch_loc = list()

    for base in reverse_trig:
        comp_base = base_pairing[base]
        base_at_loc = bstrand[counter]

        if comp_base == base_at_loc and counter + 1 <= b_base:
            pass

        elif comp_base != base_at_loc and counter + 1 <= base_count:
            base_loc_5 = b_trig - counter
            mismatch_loc.append(base_loc_5)
            mismatch += 1

        counter += 1

    overhang = b_trig - base_count
    mismatch_loc.sort()

    return mismatch, mismatch_loc, overhang


def strand_alignment(
    base_strand: str, incumbent: str, trigger: str, tb_mismatch=False, tb=7
):
    try:
        base = DNA(base_strand)
        incumb = DNA(incumbent)
        trig = DNA(trigger)
    except ValueError:
        raise TypeError("Not all sequence given is DNA sequences")

    if tb_mismatch:
        mismatch_count, mismatch_loc, ob = strand_alignment_tb_mismatch(
            base, incumb, trig, tb
        )

    else:
        mismatch_count, mismatch_loc, ob, tb = strand_alignment_no_tb_mismatch(
            base, incumb, trig
        )

    return mismatch_count, mismatch_loc, ob, tb
