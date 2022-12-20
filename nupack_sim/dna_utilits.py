import regex as re

base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}


def is_dna(seq: str):
    lower_case_seq = remove(seq).lower()
    return bool(re.match("^[atgc]+$", lower_case_seq))


def reverse_seq(seq: str):
    return seq[::-1]


def remove(seq: str):
    return seq.replace(" ", "")
