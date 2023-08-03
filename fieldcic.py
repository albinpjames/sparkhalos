"""This program is used to calculate the cic of particles
"""

# Choose the simulation to be used
from sparkhalos.simulprocess import test_rand as simulation
# from sparkhalos.simulprocess import abacussummit as simulation

from sparkhalos.simulprocess.simparams import SimuParams
from sparkhalos.hstats.cic import cic_particles
from sparkhalos.hstats.fitfun import pois, normfun


from mpl_toolkits import mplot3d
from scipy.stats import binned_statistic_dd 

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import numpy as np
import random

# The location of where the data is stored.
datalocation = "/mnt/dark/Projects/3.CUSAT/Data"

# Intilaises the simulation parameters for the simulation
params = SimuParams.init(datalocation, simulation.simparams)

# Redshifts to be computed 
# redshifts = ["3.000","2.000","1.100","0.500","0.200"]
redshifts = ["1.400"]

# The total bins taken while calculation of HMF
totalbins = 50

# The new box size & Bins for computing count in cells distribution
nw_boxsize = 25
cicbins = 15

# Choose the method for calculating the cic
cic_method = "manual"


for redshift in redshifts:
    ''' Here we are trying to calculate the CIC for various number of 
    particles generated at random.
    '''
    match params.name:
        case "test_randomnum":
            print( "Data is not read but generated on the go")
        case "abacussummit":
            field = simulation.readfieldrv(params, redshift)
            print("reading data complete")
    
    # particles_taken = [100000, 1000000, 6000000, 10000000]
    particles_taken = [100, 1000, 5000, 50000, 10000]
    cutside = int(params.boxsize/nw_boxsize)
    totalboxes = int(cutside**3) 

    # For ploting find minimium required rows given we want 2 columns 
    ncols = 2
    nrows = len(particles_taken) // ncols + 1

    # Intalizing the figure
    print("The figure is intalised")
    plt.figure(figsize=(10, 16), dpi=150)
    plt.subplots_adjust(left=0.1,
        bottom=0.1,
        right=0.9,
        top=0.8,
        wspace=0.2,
        hspace=0.5)
    plt.suptitle("CIC From Field Particles", fontsize=15, y=0.95)

    ax = plt.subplot(nrows, ncols, 1)
    ax.axis('off')
    ax.text(0.05, 0.3, 'Simulation: ' + str(params.name) +
                             '\nCosmology: LCDM' + 
                             '\nRedshift: ' + str(redshift) + 
                             '\nBox Length: ' + str(params.boxsize) + ' MPc/h' +
                             '\n\nSub Box Length: ' + str(nw_boxsize) + ' MPc/h' +
                             '\nTotal Number Of Sub Boxes: ' + str(totalboxes) +
                             '\n\nNumber Of Bins: ' + str(cicbins) +
                             '\n\nCIC method used : ' + str(cic_method),
                             fontsize=10 
                             )


    for p, n in enumerate(particles_taken): 
        # p is to find the row and column in subplot
        # n is the number of particles taken

        # Choosing the particles for which the CIC Is to be calculated
        match params.name:
            case "test_randomnum":
                data = simulation.generate_xyz(params,n)
                print("The data is simulated for test_rand")

            case "abacussummit":
                randomlist = random.sample(range(0, len(field)), n)
                data = field['pos'][randomlist]
                data += params.boxsize / 2
                print(f"Random {n} data points are choosen for abacussummit")
        
        # Choosing the CIC Calculation method
        match cic_method:
            case "manual":
                print("Calculating CIC Using Maunal Method")
                x = data[:,0:3]/nw_boxsize
                # del data
                print("converting to end")
                x = x.astype(int)

                boxes = []
                print(f"caclulating boxes {n} for particles")
                
                for i in range (len(x)):
                    
                    # Here some of the boxes maybe at the edge and hence are nudged to the adjacent boxes
                    for j in range(3):
                        if x[i][j]>= params.boxsize/nw_boxsize:
                            x[i][j] -= 1
                    
                    # This line calculates to which the boxes belongs to
                    boxes.append(x[i][0] + cutside*x[i][1] + cutside**2*x[i][2])

                #Number of particles in each boxes
                boxdata = np.zeros(totalboxes)
                for i in boxes:
                    boxdata[i] += 1
        
            case "binned_stat":
                print("Code to be updated")    
                # boxdata = binned_statistic_dd(data, statistic = 'count', bins =)
        
        # Plotting the data
        ax = plt.subplot(nrows, ncols, p+2)
        y_value, binedges, patches = ax.hist(boxdata, bins = cicbins, density=True, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha = 0.5)
        x_value = (0.5*(binedges[1:] + binedges[:-1]))
        ax.set_title(f"CIC Computed for {n} particles")
        
        # Calculating the error 
        y_valueerr, binedgeserr = np.histogram(boxdata, bins=cicbins)
        errorbar = np.sqrt(y_valueerr)/((binedges[1]-binedges[0])*np.sum(y_valueerr))
        ax.errorbar(x_value, y_value, yerr=errorbar, fmt='k.')

        # Curve Fitting Poisson
        popt_pois, pcov_pois = curve_fit(pois, x_value, y_value, p0=(np.mean(boxdata),))
        ax.plot(x_value, pois(x_value, *popt_pois), 'b-', label='Poisson (Fit): \nn=%5.3f' % tuple(popt_pois))

        # Curve Fitting Normal Dist
        popt_norm, pcov_norm = curve_fit(normfun, x_value, y_value, p0 = (np.mean(boxdata), np.std(boxdata)) )
        ax.plot(x_value, normfun(x_value, *popt_norm), 'g-', label='Normal (Fit): \nn=%5.3f, sig=%5.3f' % tuple(popt_norm))

        ax.legend(loc="upper right", fontsize=5)
        
plt.show()

def dataplot():
    """To verify the data points in space
    """
    # Creating figure
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
     

    # Creating plot
    ax.scatter3D(data[:,0],data[:,1], data[:,2], color = "green", s = 1)
    plt.title(f"For redshift {redshift}, taking {n} random points from {len(field)} data points")
     
    # show plot
    plt.show()




