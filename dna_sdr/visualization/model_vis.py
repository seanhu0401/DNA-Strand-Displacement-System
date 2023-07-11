import os, pickle
import matplotlib.pyplot as plt
import numpy as np
from dna_sdr.fitting.fit import one_phase_association, kin_fit

if __name__ == "__main__":
    os.chdir("./dna_sdr/pickles/")

    with open("kinetic_result.pkl", "rb") as fp:
        output = pickle.load(fp)
        result = output["P1_T1"]
        result.plot()
        dely = result.eval_uncertainty(sigma=3)

        # for condition in output:
        #     result = output[condition]
        #     result.plot()
        #     dely = result.eval_uncertainty(sigma=3)
