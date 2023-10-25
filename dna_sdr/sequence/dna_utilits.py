"""
Module docstring
"""

from dataclasses import dataclass, field
import regex as re


base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}


@dataclass
class DNA:
    """
    Class docstring
    """

    seq: str
    base_count: int = field(init=False)

    def __post_init__(self) -> None:
        # Check if the input sequence is DNA
        if bool(re.match("^[atgc ATGC]+$", self.seq)):
            self.seq = self.seq.replace(" ", "").lower()
            self.base_count = len(self.seq)
        else:
            raise ValueError(
                f"The sequence given is not a DNA sequence. Sequence given: {self.seq}"
            )

    def complement_seq(self) -> str:
        """
        xxx

        Parameters
        ----------

        Returns
        -------

        """
        base_list: list[str] = []
        for base in self.seq:
            base_list.append(base_pairing[base])
        complment_seq = "".join(base_list)
        return complment_seq

    def reverse_seq(self) -> str:
        """
        xxx

        Parameters
        ----------

        Returns
        -------

        """
        reverse = self.seq[::-1]
        return reverse

    def reverse_complement(self) -> str:
        """
        xxx

        Parameters
        ----------

        Returns
        -------

        """
        reversed_seq = self.reverse_seq()
        reverse_complment_seq = reversed_seq[::-1]
        return reverse_complment_seq


@dataclass
class Trigger(DNA):
    """
    Class docstring
    """

    name: str
    toehold: int
    overhang: int
    mismatch: int = 0
    mismatch_loc: list[int] | None = None
    mismatch_type: list[int] | None = None
    plate_loc: str | None = None
