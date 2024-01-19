import pytest


@pytest.fixture
def Base_Test_Seq():
    return "TCT CA CAATC CA TCT"


@pytest.fixture
def Incumb_Seq():
    return "CA CAATC CA TCT CA CCACC CA"


@pytest.fixture
def Base_Seq():
    return "TG GGTGG TG AGA TG GATTG TG AGA TG TG AGA CAT ACA GCG CCG ACC GTA"


@pytest.fixture
def Comp_Trigger():
    return "CA TAACA CA TCT CA CAATC CA TCT CA CCACC CA"
