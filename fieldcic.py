from sparkhalos.simulprocess import abacussummit
from sparkhalos.simulprocess.simparams import SimuParams
from sparkhalos.hstats.cic import cic_particles
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import random


datalocation = "/mnt/dark/Projects/3.CUSAT/Data"
params = SimuParams.init(datalocation, abacussummit.simparams, type= "small", intcont= "ph3000", boxsize=500)

# Redshifts to be computed
# redshifts = ["3.000","2.000","1.100","0.500","0.200"]
redshifts = ["1.400"]
totalbins = 50
cicbins = 10
nw_boxsize = 25

for redshift in redshifts:
    field = abacussummit.readfieldrv(params, redshift)
    print("reading data complete")
    # nw_boxsize = 25
    # print("calculating cic")
    # boxdata = cic_particles(field, params, redshift, nw_boxsize, totalbins)


    n = 1000000
    randomlist = random.sample(range(0, len(field)), n)
    data = field['pos'][randomlist]
    
    # Creating figure
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
     

    # Creating plot
    ax.scatter3D(data[:,0],data[:,1], data[:,2], color = "green", s = 1)
    plt.title(f"For redshift {redshift}, taking {n} random points from {len(field)} data points")
     
    # show plot
    plt.show()
    
    cutside = int(params.boxsize/nw_boxsize)
    totalboxes = int(cutside**3) 

    data += params.boxsize / 2
    x = data[:,0:3]/nw_boxsize
    del data
    print("converting to end")
    x = x.astype(int)

    boxes = []
    print("caclulating boxes")
    for i in range (len(x)):
        boxes.append(x[i][0] + cutside*x[i][1] + cutside**2*x[i][2])

    #Box halo numbers

    boxdata = np.zeros(totalboxes)
    for i in boxes:
        boxdata[i] += 1
        
    plt.hist(boxdata, bins = cicbins, density=True, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha = 1)
    plt.text(0.05, 0.3, 'Simulation: ' + str(params.name) +
                         '\nCosmology: LCDM' +
                         '\nRedshift: ' + str(redshift) + 
                         '\nBox Length: ' + str(params.boxsize) + ' MPc/h' +
                         '\n\nSub Box Length: ' + str(nw_boxsize) + ' MPc/h' +
                         '\nTotal Number Of Sub Boxes: ' + str(totalboxes) +
                         '\n\nNumber Of Bins: ' + str(cicbins) +
                         '\n\nNumber of Points taken: ' + str(len(randomlist)) +
                         '\n\nNumber of Points taken: ' + str(len(x))
                         )
    plt.show()

"""
    # Creating figure
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
     

    # Creating plot
    ax.scatter3D(data[:,0],data[:,1], data[:,2], color = "green", s = 1)
    plt.title(f"For redshift {redshift}, taking {n} random points from {len(field)} data points")
     
    # show plot
    plt.show()
    
import matplotlib.pyplot as plt
bins = 10
plt.hist(np.sum(boxdata, axis=0), bins=bins, density=True, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha = 0.15)

"""



