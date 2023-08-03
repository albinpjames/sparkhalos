from math import factorial as fact, sqrt as sq, exp, pi
from scipy.special import gamma, gammaln
import numpy as np


# Fitting Functions
# GQED -----------------------------------------------------------------------------------------------------
def gqed_old(x, a, b):
    return (
        a
        * (1 - b)
        / gamma(x + 1)
        * ((a * (1 - b)) + (x * b)) ** (x - 1)
        * np.exp(-a * (1 - b) - (x * b))
    )


def gqed(x, a, b):
    y = a * (1 - b) + (x * b)
    return np.exp(np.log(a * (1 - b)) - gammaln(x + 1) + (x - 1) * np.log(y) - y)


def gqed_b(n_corr):
    return 1 - (np.mean(n_corr) + 1) ** (-0.5)


# Poisson --------------------------------------------------------------------------------------------------
def pois(x, a):
    return np.exp(x * np.log(a) - a - gammaln(x + 1))


# Gaussian -------------------------------------------------------------------------------------------------
def normfun(x, sig, a):
    return (1 / (sig * sq(2 * pi))) * np.exp(-0.5 * ((x - a) / sig) ** 2)




