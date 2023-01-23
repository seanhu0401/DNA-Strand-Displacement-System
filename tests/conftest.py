import pytest
from dna_sdr.sequence.dna_utilits import DNA


@pytest.fixture
def Base_Test_Seq():
    return "TCT CA CAATC CA TCT"


@pytest.fixture
def Incumb_Seq():
    return "CA CAATC CA TCT CA CCACC CA"


@pytest.fixture
def Base_Seq():
    return "TG GGTGG TG AGA TG GATTG TG AGA TG TG AGA CAT ACA GCG CCG ACC GTA"
