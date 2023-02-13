import numpy as np


def linear_function(x, slope, intercept):
    return slope * x + intercept


def logistic_function(x, L, k, xo):
    return L / (1 + np.exp(-k * (x - xo)))
