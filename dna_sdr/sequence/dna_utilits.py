import regex as re
from dataclasses import dataclass, field
from typing import Optional

base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}


@dataclass
class DNA:
    seq: str
    base_count: int = field(init=False)

    def __post_init__(self):
        # Check if the input sequence is DNA
        if bool(re.match("^[atgc ATGC]+$", self.seq)):
            self.seq = self.seq.replace(" ", "").lower()
            self.base_count = len(self.seq)
        else:
            raise ValueError(
                "The sequence given is not a DNA sequence. Sequence given: {}".format(
                    self.seq
                )
            )

    def complement(self):
        base_list = list()
        for base in self:
            base_list.appen(base_pairing[base])
        complment_seq = "".join(base_list)
        return complment_seq

    def reverse_seq(self):
        reverse = self.seq[::-1]
        return reverse

    def reverse_complement(self):
        reversed = self.reverse_seq()
        reverse_complment_seq = self.complement(reversed)
        return reverse_complment_seq


@dataclass
class trigger(DNA):
    name: str
    toehold: Optional[int] = None
    overhang: Optional[int] = None
    mismatch: Optional[int] = None
    mismatch_loc: Optional[list] = None
    mismatch_type: Optional[list] = None
    plate_loc: Optional[str] = None

    def __post_init__(self):
        super().__post_init__()
