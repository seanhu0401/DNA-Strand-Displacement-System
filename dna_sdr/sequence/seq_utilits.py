from dna_sdr.sequence.dna_utilits import DNA

comp_dict = {"t": 1, "g": 2, "c": 3, "a": 4}
base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}

# *Think of better name for the method
# This function is currently used to compare to sequence, mainly the triggers, to identify the location(s) that has/have mismatch(es)
# They should be equal length in bases if there is no modification to the toehold + dangles - not applicable with those modification currently
def sequence_comparison(seq_1: DNA, seq_2: DNA):
    comparison_lst = []
    base_loc_lst = []
    counter = 0

    if type(seq_1) != DNA and type(seq_2) != DNA:
        raise TypeError("The input is not a DNA type")

    elif seq_1.base_count != seq_2.base_count:
        raise Exception("The two sequences are not equal length")

    else:
        while counter < seq_1.base_count:
            comp_val = 0
            if seq_1.seq[counter] == seq_2.seq[counter]:
                comparison_lst.append(comp_val)
            elif seq_1.seq[counter] != seq_2.seq[counter]:
                comp_val = comp_dict[seq_2.seq[counter]]
                comparison_lst.append(comp_val)
            base_loc_lst.append(counter + 1)
            counter += 1
        seq_comp_dict = dict(zip(base_loc_lst, comparison_lst))

    return seq_comp_dict


def strand_alignment(
    base_strand: DNA,
    incumbent: DNA,
    trigger: DNA,
    toehold_base_number=None,
    toehold_base_mismatch=False,
):
    if type(base_strand) != DNA:
        raise TypeError(
            "The sequence given is not a DNA sequence. Sequence given: {}".format(
                base_strand
            )
        )

    elif type(incumbent) != DNA:
        raise TypeError(
            "The sequence given is not a DNA sequence. Sequence given: {}".format(
                incumbent
            )
        )
    elif type(trigger) != DNA:
        raise TypeError(
            "The sequence given is not a DNA sequence. Sequence given: {}".format(
                trigger
            )
        )

    bstrand = base_strand.seq()
    b_base = base_strand.base_count()

    incumb = incumbent.seq()
    b_incumb = incumbent.base_count()

    reverse_trig = trigger.reverse_seq()
    b_trig = trigger.base_count()

    if not toehold_base_mismatch:
        counter = total_mismatch_count = last_match_base_loc = 0
        mismatch_loc = list()
        tb_ob = b_trig - b_incumb  # number of bases for toehold + overhang

        for base in reverse_trig:
            comp_base = base_pairing[base]
            base_at_loc = bstrand[counter]
            if comp_base == base_at_loc and counter <= b_base:
                pass

            elif comp_base != base_at_loc and b_trig - counter >= tb_ob:
                base_loc_5 = b_trig - counter  # Base location from the 5' end
                mismatch_loc.append(base_loc_5)
                total_mismatch_count += 1

            elif (
                comp_base != base_at_loc
                and b_trig - counter <= tb_ob
                and counter > last_match_base_loc
            ):
                if last_match_base_loc < b_trig - counter:
                    last_match_base_loc = counter
            counter += 1

        if last_match_base_loc != 0:
            toehold = last_match_base_loc - b_incumb
            overhang = tb_ob - toehold

        else:
            toehold = tb_ob
            overhang = 0

    else:
        counter = total_mismatch_count = 0
        mismatch_loc = list()
        total_base_number = b_incumb + toehold_base_number

        for base in reverse_trig:
            comp_base = base_pairing[base]
            base_at_loc = bstrand[counter]

            if comp_base == base_at_loc and counter + 1 <= (total_base_number):
                pass

            elif counter + 1 > (total_base_number):
                overhang = b_trig - counter
                break

            else:
                base_loc_5 = b_trig - counter
                mismatch_loc.append(base_loc_5)
                total_mismatch_count += 1

            counter += 1

        toehold = toehold_base_number

    return mismatch_loc, overhang, toehold
