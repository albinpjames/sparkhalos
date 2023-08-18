from sparkhalos.hstats.fitfun import gev, gevold

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import numpy as np

x_value = np.arange(-1,3,0.1)
y_value = gev(x_value,1,1,1)
plt.plot(x_value, y_value)
plt.show()
