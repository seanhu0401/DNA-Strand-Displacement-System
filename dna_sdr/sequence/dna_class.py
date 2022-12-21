from dataclasses import dataclass, field
import dna_utilits as utilis

base_pairing = {"a": "t", "t": "a", "g": "c", "c": "g"}


@dataclass
class DNA:
    __slot__ = ["seq", "base_count"]
    seq: str
    base_count: int = field(init=False)

    def __post_init__(self):
        self.seq = utilis.remove(self.seq).lower()
        self.base_count = len(self.seq)

    def reverse_complement(self):
        inverse = utilis.reverse_seq(self.seq)
        base_list = list()
        for base in inverse:
            base_list.append(base_pairing[base])
        complment_seq = "".join(base_list)
        return complment_seq


# TODO: Make a dataclass for triggers
# @dataclass
# class trigger(DNA):


if __name__ == "__main__":
    base = DNA("TG GGTGG TG AGA TG GATTG TG AGA TG TG AGA CAT ACA GCG CCG ACC GTA")
    print(base.base_count)
    print(base.seq[1])
