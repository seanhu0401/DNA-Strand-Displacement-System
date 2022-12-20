import numpy as np


def one_phase_association(time, plateau, k):
    return plateau * (1 - np.exp(-k * time))
