import pandas as pd
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def simple_reaction_kinetic(X, t, const):
    k1, k2 = const
    I = X[0]
    QT = X[1]
    IT = X[2]
    Q = X[3]

    dIdt = -k1 * I * QT + k2 * IT * Q
    dQTdt = -k1 * I * QT + k2 * IT * Q
    dITdt = k1 * I * QT - k2 * IT * Q
    dQdt = k1 * I * QT - k2 * IT * Q

    f = np.array([(dIdt), (dQTdt), (dITdt), (dQdt)])

    return f


def complex_reaction_kinetic(X, t, const):
    k1, k2 = const
    I = X[0]
    QT = X[1]
    IQT = X[2]
    IT = X[3]
    Q = X[4]

    dIdt = -k1 * IQT * I
    dQTdt = -k1 * IQT * I
    dIQTdt = k1 * I * QT - k2 * IQT
    dITdt = k2 * IQT
    dQdt = k2 * IQT

    f = np.array([(dIdt), (dQTdt), (dIQTdt), (dITdt), (dQdt)])

    return f


pickle = "dna_sdr/IO/Output/Pickles/4WJ_HEX_Screen_P0_A01_A04_summerized.pkl"
df = pd.read_pickle(pickle)

t1 = df["T1_mean"]
time = df["time (min)"]

C0 = np.zeros(4, dtype=float)
C0[0] = 500
C0[1] = 500

k1 = 5e-3  # nM/min
k2 = 0
const = [k1, k2]

soln = odeint(simple_reaction_kinetic, C0, time, args=(const,))
I, QT, IT, Q = soln.T

print((t1 * 500 - IT) ** 2)

print(sum((t1 * 500 - IT) ** 2))

# soln = odeint(complex_reaction_kinetic, C0, time, args=(const,))
# I, QT, IQT, IT, Q = soln.T

# plt.plot(time, I, color="red", label="Trig")
# plt.plot(time, QT, color="green", label="Inter")
# plt.scatter(time, t1 * 500, label="Exp")
# plt.plot(time, IT, color="blue", label="Product")
# plt.legend()

# plt.tight_layout()
# plt.show()
