import regex as re
from dataclasses import dataclass, field

base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}


@dataclass
class DNA:
    __slot__ = ["seq", "base_count"]
    seq: str
    base_count: int = field(init=False)

    def __post_init__(self):
        # *Check if the input sequence is DNA
        if bool(re.match("^[atgc ATGC]+$", self.seq)):
            self.seq = self.seq.replace(" ", "").lower()
            self.base_count = len(self.seq)
        else:
            raise ValueError(
                "The sequence given is not a DNA sequence. Sequence given: {}".format(
                    self.seq
                )
            )

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
@dataclass
class trigger(DNA):
    toehold: int
