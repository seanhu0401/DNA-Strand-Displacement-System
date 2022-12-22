import pytest
from dna_sdr.sequence.dna_utilits import DNA


@pytest.fixture
def Base_Test_Seq():
    return DNA("TCT CA CAATC CA TCT")
