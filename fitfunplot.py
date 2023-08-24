"""This program is used to calculate the cic of particles
"""

# Choose the simulation to be used
# from sparkhalos.simulprocess import test_rand as simulation
from sparkhalos.simulprocess import abacussummit as simulation

from sparkhalos.simulprocess.simparams import SimuParams
from sparkhalos.hstats.cic import cic_particles
from sparkhalos.hstats.fitfun import pois, normfun, gev, lnnorm

from scipy.stats import binned_statistic_dd 
from scipy.stats import skew
from scipy.stats import genextreme


import matplotlib.pyplot as plt
import scienceplots
plt.style.use('science')
from scipy.optimize import curve_fit

import numpy as np


# The location of where the data is stored.
datalocation = "/mnt/dark/Projects/3.CUSAT/Data/3cic2000"
# datalocation = "/home/darkmatter/Documents/Albin/DATA/3cic2000"
# datalocation = "/home/albinpjames/Documents/CUSAT/DATA/3cic2000"

params = SimuParams.init(datalocation, simulation.simparams)

nw_boxsizes = [10,20,25,50]
cicbins = 20
size = 2000
redshift = 3

ncols = 2
nrows = len(nw_boxsizes) // ncols + 1

# Intalizing the figure
print("The figure is intalised")
plt.figure(figsize=(8, nrows * 3), dpi=150)
plt.subplots_adjust(left=0.1,
    bottom=0.1,
    right=0.9,
    top=0.8,
    wspace=0.2,
    hspace=0.5)
plt.suptitle("CIC From Field Particles", fontsize=15, y=0.95)

# ax = plt.subplot(nrows, ncols, 1)
# ax.axis('off')
# ax.text(0.05, 0.3, 'Simulation: ' + str(params.name) +
#                          '\nCosmology: LCDM' + 
#                          '\nRedshift: ' + str(redshift) + 
#                          '\nBox Length: ' + str(params.boxsize) + ' MPc/h' +
#                          # '\n\nSub Box Length: ' + str(nw_boxsize) + ' MPc/h' +
#                          # '\nTotal Number Of Sub Boxes: ' + str(totalboxes) +
#                          # '\nNumber of particles taken: ' + str(particles_taken) +
#                          '\n\nNumber Of Bins: ' + str(cicbins) +
#                          # '\n\nCIC method used : ' + str(cic_method),
#                          fontsize=10 
#                          )

for p, nw_boxsize in enumerate(nw_boxsizes): 
	
	particles_in_box = np.load(datalocation + "/" + str(nw_boxsize) +".npy")
	print(f"Calculating for boxsize {nw_boxsize} with particles: {np.sum(particles_in_box)}")

	
	'Here we are converting the particles in box to density contrast '
	cell_avg = np.sum(particles_in_box) * (nw_boxsize**3) / (size**3)
	celldensity = (particles_in_box - cell_avg)/cell_avg
	boxdata = celldensity # The density contrast
	boxdata = np.log(boxdata+1)

	'Plotting the data'
	ax = plt.subplot(nrows, ncols, p+1)
	ax.set_title(f"CIC Computed for {nw_boxsize} sub box size")
	# ax.set_xscale('log')

	'The data is binned and x and y values are obtained'
	hist_y, hist_edge = np.histogram(boxdata, bins=cicbins)

	normfactor = ((hist_edge[1]-hist_edge[0])*np.sum(hist_y))
	y_value = hist_y / normfactor
	
	'The delta is converted to log(delta+1)'
	x_center = (0.5*(hist_edge[1:] + hist_edge[:-1]))
	# x_value = np.log(x_center+1)
	x_value = x_center

	# ax.scatter(x_value,y_value)

	'Error'
	errorbar = np.sqrt(hist_y) / normfactor
	ax.errorbar(x_value, y_value, yerr=errorbar, fmt='k.')

		# Curve Fitting GEV
	print(f"Calculating fitting for GEV {nw_boxsize}")
	popt_gev, pcov_gev = curve_fit(gev, x_value, y_value, 
		# p0 =(skew(boxdata),np.mean(boxdata),np.std(boxdata)),
		p0 =(skew(boxdata),np.std(boxdata)),
		)
	print(f"Plotting GEV {nw_boxsize}")
	# xvalue = np.arange(-1,1,0.1 )
	xvalue = x_value
	# ax.plot(xvalue, gev(xvalue, *popt_gev), 'g--', label='GEV (Fit): \nxi=%5.3f, \nnu_g=%5.3f, \nsig_g=%5.3f' % tuple(popt_gev))
	ax.plot(xvalue, gev(xvalue, *popt_gev), 'g--', label='GEV (Fit): \nxi=%5.3f, \nnu_g= 0, \nsig_g=%5.3f' % tuple(popt_gev))

	# # Curve Fitting GEV SciP
	# print(f"Calculating fitting for Genextreme {nw_boxsize}")
	# # popt_genex = genextreme.fit(x_value)
	# popt_genex, pcov_gev = curve_fit(genextreme.pdf, x_value, y_value, p0 =(skew(boxdata),np.mean(boxdata),np.std(boxdata)))
	# xi, nu_g, sig_g = popt_genex
	# print(f"Plotting Genextreme {nw_boxsize}")
	# ax.plot(x_value, genextreme.pdf(x_value, - xi, nu_g, sig_g), 'b--', label='Genextreme (Fit): \nxi=%5.3f, \nnu_g=%5.3f, \nsig_g=%5.3f' % tuple(popt_genex))


    # Curve Fitting ln Norm
	# print(f"Calculating fitting for Log Normal {nw_boxsize}")
	# popt_lnnrm, pcov_lnnrm = curve_fit(lnnorm, x_value, y_value, p0 =(np.mean(boxdata),np.std(boxdata)))
	# print(f"Plotting Log Normal {nw_boxsize}")
	# ax.plot(x_value, lnnorm(x_value, *popt_lnnrm), 'r--', label='Lognorm (Fit): \nnu =%5.3f, \nsig =%5.3f' % tuple(popt_lnnrm))
	
	# xvalue = np.arange(-1,4,0.1 )
	# ax.plot(xvalue, genextreme.pdf(xvalue,-0.5,0,0.4), label='scipy genextreme')
	# ax.plot(xvalue, gev(xvalue, 0.5,0,0.4), label='gev')


	# Curve Fitting Poisson
	# popt_pois, pcov_pois = curve_fit(pois, x_value, y_value, p0=(np.mean(boxdata),))
	# ax.plot(x_value, pois(x_value, np.mean(boxdata)), 'b--', label=f'Mean (Poisson Fit):{popt_pois} \n Actual Mean {np.mean(boxdata)}' )

	# Curve Fitting Normal Dist
	# popt_norm, pcov_norm = curve_fit(normfun, x_value, y_value, p0 = (np.mean(boxdata), np.std(boxdata)) )
	# ax.plot(x_value, normfun(x_value, *popt_norm), 'b--', label='Normal (Fit): \nn=%5.3f, sig=%5.3f' % tuple(popt_norm))

	ax.legend(loc="upper right", fontsize=4)

plt.show()