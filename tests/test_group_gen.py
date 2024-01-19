import pytest
import pandas as pd
from dna_sdr.experimental.data_processing import group_generation


def conversion(number_list):
    str_list = list()
    for i in number_list:
        s = list(map(str, i))
        str_list.append(s)
    return str_list


test_lst = [
    (
        "./dna_sdr/IO/Output/Pickles/4WJ_HEX_Screen_P0_A01_A04_normalized.pkl",
        21,
        [
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12],
            [13, 14, 15, 16],
            [17, 18, 19, 20],
        ],
        [
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12],
            [13, 14, 15, 16],
            [17, 18, 19, 20],
        ],
    ),
    (
        "./dna_sdr/IO/Output/Pickles/4WJ_HEX_Screen_P0_A05_A08_normalized.pkl",
        21,
        [
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12],
            [13, 14, 15, 16],
            [17, 18, 19, 20],
        ],
        [
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12],
            [13, 14, 15, 16],
            [17, 18, 19, 20],
        ],
    ),
    (
        "./dna_sdr/IO/Output/Pickles/4WJ_HEX_Screen_P1_A11_B03_normalized.pkl",
        25,
        [
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12],
            [13, 14, 15, 16],
            [17, 18, 19, 20],
            [21, 22, 23, 24],
        ],
        [
            [1, 2, 3, 4],
            [5, 6, 7, 8],
            [9, 10, 11, 12],
            [13, 14, 15, 16],
            [17, 18, 19, 20],
            [21, 22, 23, 24],
        ],
    ),
]


@pytest.mark.parametrize("input, control_number, sample, presented", test_lst)
def test_group_gen(input, control_number, sample, presented):
    df = pd.read_pickle(input)
    output = group_generation(df, control_number)
    str_sample = conversion(sample)
    str_presented = conversion(presented)

    assert str_sample == output[0]
    assert str_presented == output[1]


if __name__ == "__main__":
    df = pd.read_pickle(
        "./dna_sdr/IO/Output/Pickles/4WJ_HEX_Screen_P0_A01_A04_normalized.pkl"
    )
    print(df)
