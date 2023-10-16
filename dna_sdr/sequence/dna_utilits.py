from dataclasses import dataclass, field
from typing import Optional
import regex as re


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
                f"The sequence given is not a DNA sequence. Sequence given: {self.seq}"
            )

    def complement(self):
        base_list = []
        for base in self.seq:
            base_list.append(base_pairing[base])
        complment_seq = "".join(base_list)
        return complment_seq

    def reverse_seq(self):
        reverse = self.seq[::-1]
        return reverse

    def reverse_complement(self):
        reversed_seq = self.reverse_seq()
        reverse_complment_seq = reversed_seq.complement()
        return reverse_complment_seq


@dataclass
class Trigger(DNA):
    name: str
    toehold: Optional[int] = None
    overhang: Optional[int] = None
    mismatch: Optional[int] = None
    mismatch_loc: Optional[list] = None
    mismatch_type: Optional[list] = None
    plate_loc: Optional[str] = None
