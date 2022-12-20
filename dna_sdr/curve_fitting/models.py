import numpy as np

"""
The equation(s) below was used to fit the release curve from the qPCR experiments.
"""

# pseudo-first order association kinetics
def one_phase_association(time, plateau, k):
    return plateau * (1 - np.exp(-k * time))
