import regex as re
from dataclasses import dataclass, field

base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}


def remove(seq: str):
    return seq.replace(" ", "")


def is_dna(seq: str):
    lower_case_seq = remove(seq).lower()
    return bool(re.match("^[atgc]+$", lower_case_seq))


def reverse_seq(seq: str):
    return seq[::-1]


@dataclass
class DNA:
    __slot__ = ["seq", "base_count"]
    seq: str
    base_count: int = field(init=False)

    def __post_init__(self):
        self.seq = remove(self.seq).lower()
        self.base_count = len(self.seq)

    def reverse_seq(self):
        self.reverse = self.seq[::-1]

    def reverse_complement(self):
        inverse = self.reverse
        base_list = list()
        for base in inverse:
            base_list.append(base_pairing[base])
        complment_seq = "".join(base_list)
        return complment_seq


# TODO: Make a dataclass for triggers
# @dataclass
# class trigger(DNA):
