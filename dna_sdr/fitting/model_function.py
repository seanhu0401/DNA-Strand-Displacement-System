"""_summary_

Returns
-------
_type_
    _description_
"""

import numpy as np
from lmfit.lineshapes import logistic
from scipy.integrate import solve_ivp


def one_phase_association(time: float, plateau: float, k: float, y_0: float) -> float:
    """The function calculates the value of a one-phase association reaction over time."""
    return y_0 + (plateau - y_0) * (1 - np.exp(-k * time))


def lag_one_phase_association(
    time: float, plateau: float, k: float, time_0: float, y_0: float
) -> float:
    """The function calculates the value of a one-phase association reaction over time.

    Parameters
    ----------
    time : float
        The time parameter represents the time at which the association is being calculated.
    plateau : float
        The plateau is the maximum value that the function will approach as time goes to infinity.
    k : float
        The parameter "k" represents the rate constant in the one-phase association equation. It determines
    how quickly the association reaction occurs. A higher value of "k" indicates a faster association
    rate.
    time_0 : float
        The value where the reaction start after the initial lag.
    y_0 : float
        The initial value during the initial lag time.

    Returns
    -------
        y_0 if the time is less than time_0 else result of the equation
        `y_0 + (plateau - y_0) * (1 - np.exp(-k * time))`.

    """
    step = logistic(time, center=time_0, sigma=0.1)
    return (
        y_0 * (1 - step)
        + (y_0 + (plateau - y_0) * (1 - np.exp(-k * (time - time_0)))) * step
    )


def second_kinetic(t: list[float], y0: list[float], k: float) -> list[float]:
    """The function `second_kinetic` defines a simplified system of ordinary differential equations (ODEs)
    for a second order irreversible reaction.

    Parameters
    ----------
    t : list[float]
        The parameter `t` is a list of time values at which you want to evaluate the ODEs. It represents
    the time points at which you want to calculate the values of the variables in the system.
    y0 : list[float]
        The parameter `y0` is a list of initial values for the variables in the system. In this case, the
    variables are `I`, `QT`, `IT`, and `Q`. So `y0` should be a list of four floats representing the
    initial values of these variables.
    k : float
        The parameter `k` represents the rate constant for the second order irreversible reaction. It
    determines the rate at which the reaction occurs.

    Returns
    -------
        a list of floats, which represents the system of ordinary differential equations (ODEs) for the
    given second order irreversible reaction.

    """

    I, QT, IT, Q = y0
    k1 = k

    # the model equations
    dIdt = -k1 * I * QT
    dQTdt = -k1 * I * QT
    dITdt = k1 * I * QT
    dQdt = k1 * I * QT
    ode = [dIdt, dQTdt, dITdt, dQdt]
    return ode


def second_kinetic_IVP_solver(t: list[float], y0: list[float], k1: float) -> np.ndarray:
    """The function `second_kinetic_IVP_solver` solves a second-order kinetic equation using the `solve_ivp` function and
    returns the third element of the solution.

    Parameters
    ----------
    t : list[float]
        A list of time values at which the solution is evaluated.
    y0 : list[float]
        The parameter `y0` represents the initial conditions for the system of differential equations. In
    this case, it is a list of initial values for the dependent variables in the system. The length of
    `y0` should match the number of equations in the system.
    k1 : float
        The parameter `k1` represents the rate constant for the second-order kinetic reaction. It
    determines the rate at which the reaction proceeds.

    Returns
    -------
        The function `sec_kin_fit` returns an `np.ndarray` which represents the third element of the
    solution of the second order kinetic equation.

    """
    release_matrix = solve_ivp(
        second_kinetic,
        (0, max(t) + 10),
        y0,
        t_eval=t,
        args=(k1,),
        rtol=1e-8,
        method="RK45",
    )
    return release_matrix.y[2]
