"""This program is used to calculate the cic of particles
"""

import sys
sys.path.append('..')

from sparkhalos.simparams.simparams import SimuParams
""" Choose the simulation to be used """
# from sparkhalos.simulations import test_rand as simulation
# from sparkhalos.simparams.test_rand import test500 as simparams
from sparkhalos.simulations import abacussummit as simulation
from sparkhalos.simparams.abacussummit import hugebase2000 as simparams

from localfiles import getlocation, savelocation


from sparkhalos.hstats.fitfun import pois, normfun, gev_delta, lnnorm, gev_A

from scipy.stats import binned_statistic_dd 
from scipy.stats import skew
from scipy.stats import genextreme


import matplotlib.pyplot as plt
import scienceplots
plt.style.use('science')
from scipy.optimize import curve_fit

import numpy as np
import pandas as pd

if __name__ == "__main__":

	redshift = "3.000"

	""" Intilaises the simulation parameters for the simulation """
	datalocation = getlocation()
	params = SimuParams.init(datalocation, simparams)
	datalocation = savelocation(params,redshift)


	nw_boxsizes = [25,15,10,5,3,1]
	cicbins = 20
	size = params.boxsize
	

	ncols = 2
	nrows = len(nw_boxsizes) // ncols + 1

	dmin = 0.0
	dmax = 0.2
	# Intalizing the figure
	print("The figure is intalised")
	plt.figure(figsize=(8, nrows * 3), dpi=150)
	plt.subplots_adjust(left=0.1,
	    bottom=0.1,
	    right=0.9,
	    top=0.8,
	    wspace=0.2,
	    hspace=0.5)
	plt.suptitle(f"Halos CIC of {params.name} of size {params.boxsize} MPc/h in delta range {dmin} to {dmax}", fontsize=15, y=0.95)

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

    boxbins = np.linspace(0,params.boxsize,5)
    boxbins = boxbins[0:-1]
    dbin = boxbins[1] - boxbins[0]

	for p, nw_boxsize in enumerate(nw_boxsizes): 
		nwboxloc = os.path.join(pathsave,str(nw_boxsize))
		for pos_start in itertools.product(boxbins,boxbins,boxbins):
			cicdata = np.load(datalocation + "/" + str(nw_boxsize)  + "/" + str(pos_start) +".npy")
			massdata = np.load(datalocation + "/" + str(nw_boxsize) + "/"  + "mass_bins" +".npy")

			colm = ["M{}".format(i) for i in range(1,21)]
			massbins = ["M{}".format(i) for i in range(1,21)]
			colm.append("P")
			colm.append("D")
			boxdata = pd.DataFrame(cicdata, columns = colm)

			pardata = sum(boxdata["P"])
			print(f"Calculating for boxsize {nw_boxsize} with particles: {pardata:.1e}")

		
		'Plotting the data'
		ax = plt.subplot(nrows, ncols, p+1)
		ax.set_title(f"Sub box size: {nw_boxsize}")
		# ax.set_xscale('log')



		print(boxdata[massbins][(boxdata['D']>dmin) & (boxdata['D']<dmax)])
		testboxdata =  boxdata
		boxdata = boxdata[massbins][(boxdata['D']>dmin) & (boxdata['D']<dmax)].sum(axis=1)
		print(boxdata)
		# boxdata = boxdata["D"]
		
		'The data is binned and x and y values are obtained'
		mindata = min(boxdata)
		maxdata = max(boxdata)
		if (maxdata - mindata) <=20:
			step = 1
		else:	
			step =int((maxdata - mindata)/20)
		cicbins = np.arange(mindata -0.5, maxdata +0.5, step)
		# cicbins = 30

		hist_y, hist_edge = np.histogram(boxdata, bins=cicbins)

		normfactor = ((hist_edge[1]-hist_edge[0])*np.sum(hist_y))

		print(f"hist sum {np.sum(hist_y)} - normfactor: {normfactor}")
		y_value = hist_y / normfactor
		
		x_center = (0.5*(hist_edge[1:] + hist_edge[:-1]))
		x_value = x_center

		# ax.scatter(x_value,y_value)

		'Error'
		errorbar = np.sqrt(hist_y) / normfactor
		ax.errorbar(x_value, y_value, yerr=errorbar, fmt='k.')

		 # Curve Fitting GEV
		# print(f"Calculating fitting for GEV {nw_boxsize}")
		# popt_gev, pcov_gev = curve_fit(gev_A, x_value, y_value, 
		# 	p0 =(skew(boxdata),np.mean(boxdata),np.std(boxdata)))
		
		# print(f"Plotting GEV {nw_boxsize}")
		# xvalue = np.arange(-1,1,0.1 )
		# xvalue = x_value
		# ax.plot(xvalue, gev_A(xvalue, *popt_gev), 'g--', label='GEV (Fit): \nxi=%5.3f, \nnu_g=%5.3f, \nsig_g=%5.3f' % tuple(popt_gev))
		# ax.plot(xvalue, gev_delta(xvalue, *popt_gev), 'g--', label='GEV (Fit): \nxi=%5.3f, \nnu_g=%5.3f, \nsig_g=%5.3f' % tuple(popt_gev))

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
		popt_pois, pcov_pois = curve_fit(pois, x_value, y_value, p0=(np.mean(boxdata)))
		ax.plot(x_value, pois(x_value, np.mean(boxdata)), 'b--', label=f'Mean (Poisson Fit):{popt_pois} \n Actual Mean {np.mean(boxdata)}' )

		# Curve Fitting Normal Dist
		# popt_norm, pcov_norm = curve_fit(normfun, x_value, y_value, p0 = (np.mean(boxdata), np.std(boxdata)) )
		# ax.plot(x_value, normfun(x_value, *popt_norm), 'b--', label='Normal (Fit): \nn=%5.3f, sig=%5.3f' % tuple(popt_norm))

		ax.legend(loc="upper right", fontsize=4)

	plt.show()