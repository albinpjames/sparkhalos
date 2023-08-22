"""This program is used to calculate the cic of particles
"""

# Choose the simulation to be used
# from sparkhalos.simulprocess import test_rand as simulation
from sparkhalos.simulprocess import abacussummit as simulation

from sparkhalos.simulprocess.simparams import SimuParams
from sparkhalos.hstats.cic import cic_particles
from sparkhalos.hstats.fitfun import pois, normfun, gev, gev_mod

from scipy.stats import binned_statistic_dd 
from scipy.stats import skew
from scipy.stats import genextreme


import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import numpy as np


# # The location of where the data is stored.
# # datalocation = "/mnt/dark/Projects/3.CUSAT/Data"
# # datalocation = "/home/darkmatter/Documents/Albin/DATA"
# datalocation = "/home/albinpjames/Documents/CUSAT/DATA/3cic2000"

# nw_boxsize = 20
# cicbins = 20
# boxdata = np.load(datalocation + "/" + str(nw_boxsize) +".npy")

# size = 2000
# cell_avg = np.sum(boxdata) * (nw_boxsize**3) / (size**3)
# celldensity = (boxdata - cell_avg)/cell_avg
# oldboxdata = boxdata
# boxdata = celldensity

# # Plotting the data
ax = plt

# y_value, binedges, patches = ax.hist(boxdata, bins = cicbins, density=True, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha = 0.5)
# x_value = (0.5*(binedges[1:] + binedges[:-1]))
# ax.suptitle(f"CIC Computed for {nw_boxsize} sub box size")
# # x_value = abs(x_value)

# # Calculating the error 
# y_valueerr, binedgeserr = np.histogram(boxdata, bins=cicbins)
# errorbar = np.sqrt(y_valueerr)/((binedges[1]-binedges[0])*np.sum(y_valueerr))
# ax.errorbar(x_value, y_value, yerr=errorbar, fmt='k.')

# # Curve Fitting GEV SciPy
# popt_genex = genextreme.fit(x_value)
# print(popt_genex)
# ax.plot(x_value, genextreme.pdf(x_value, *popt_genex), 'b--', label='genextreme (Fit): xi=%5.3f, \nnu_g=%5.3f, sig_g=%5.3f' % tuple(popt_genex))


# # Curve Fitting GEV
# popt_gev, pcov_gev = curve_fit(gev, x_value, y_value, p0 =(np.mean(boxdata),np.std(boxdata),skew(boxdata,bias=True)))
# # popt_gev, pcov_gev = curve_fit(gev, x_value, y_value)
# ax.plot(x_value, gev(x_value, *popt_gev), 'g--', label='GEV (Fit): \nnu_g=%5.3f, sig_g=%5.3f, xi=%5.3f' % tuple(popt_gev))

xvalue = np.arange(-1,4,0.1 )
ax.plot(xvalue, genextreme.pdf(xvalue,-0.5,0,0.4), label='scipy genextreme')
ax.plot(xvalue, gev(xvalue, 0,0.4,0.5), label='gev')

# Curve Fitting Poisson
# popt_pois, pcov_pois = curve_fit(pois, x_value, y_value, p0=(np.mean(boxdata),))
# ax.plot(x_value, pois(x_value, np.mean(boxdata)), 'b--', label=f'Mean (Poisson Fit):{popt_pois} \n Actual Mean {np.mean(boxdata)}' )

# Curve Fitting Normal Dist
# popt_norm, pcov_norm = curve_fit(normfun, x_value, y_value, p0 = (np.mean(boxdata), np.std(boxdata)) )
# ax.plot(x_value, normfun(x_value, *popt_norm), 'g-', label='Normal (Fit): \nn=%5.3f, sig=%5.3f' % tuple(popt_norm))

ax.legend(loc="upper right", fontsize=5)

plt.show()