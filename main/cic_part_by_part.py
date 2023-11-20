"""This program is used to calculate the CIC and various parameters
"""

import sys
sys.path.append('..')

from sparkhalos.simparams.simparams import SimuParams

""" Choose the simulation to be used """
from sparkhalos.simulations import test_rand as simulation
from sparkhalos.simparams.test_rand import test100 as simparams

# from sparkhalos.simulations import abacussummit as simulation
# from sparkhalos.simparams.abacussummit import hugebase2000 as simparams

from sparkhalos.hstats.cic import cic, dens_contrast
from localfiles import getlocation, savelocation
import numpy as np
import math
import os
import itertools

def dataplot(data):
    """To verify the data points in space
    """
    from mpl_toolkits import mplot3d
    # Creating figure
    fig = plt.figure(data3dplot, figsize = (10, 7))
    ax = plt.axes(projection ="3d")
     

    # Creating plot
    ax.scatter3D(data[:,0],data[:,1], data[:,2], color = "green", s = 1)
    plt.title("Scatter Plot to test data points")
     
    # show plot
    print("showing plot")
    plt.show()

if __name__ == "__main__":

    """ Intilaises the simulation parameters for the simulation """
    datalocation = getlocation()
    params = SimuParams.init(datalocation, simparams)

    """ Redshifts & boxsizes to be computed """
    # redshifts = ["3.000","2.500","2.000"]
    redshifts = ["3.000"]
    # nw_boxsizes = [30]
    nw_boxsizes = [25,15,10,5,3,1]


    # cic_method = "manual"
    cic_method = "binned_stat"
    halobins = 20
     
    """ Number of particles taken or generated """
    particles_taken = 10**8

    halos_taken = 10**6
    take_all = True

    for redshift in redshifts:
        print(f"redshift being computed: {redshift}")
        """ Here the data is read to be processed """
        match params.name:
            case "test_randomnum":
                halodata = simulation.generate_mxyz(params,halos_taken)
                halopartno = halodata["N"]
                halodata["N"] = halodata["N"] * params.mass
                partdata = simulation.generate_xyz(params,particles_taken)
                print("The data is simulated for test_rand")

            case "abacussummit":
                """Read halo data"""
                print(f"Halo data is being read for redshift:{redshift}")
                halodata = simulation.mass_pos(params, redshift, mode="all")
                halopartno = halodata["N"]
                halodata["N"] = halodata["N"] * params.mass


                if halos_taken <= len(halodata) or take_all == False:
                    randomlist = np.random.choice(len(halodata), halos_taken, replace=False)
                    halodata = halodata['N','xpos','ypos','zpos'][randomlist]

                """Read particle data"""
                print(f"Particle data is being read for redshift:{redshift}")
                partdata = simulation.read_particles(params, redshift, [["field","halo"], ["A"]])
                print(f"No of particles in redshift- {redshift}: {len(partdata)}.")

                if particles_taken <= len(partdata) or take_all == False:
                    randomlist = np.random.choice(len(partdata), particles_taken, replace=False)
                    partdata = partdata['xpos','ypos','zpos'][randomlist]


        # The binedges for dividing the larger box to smaller boxes
        boxbins = np.linspace(0,params.boxsize,5)
        boxbins = boxbins[0:-1]
        dbin = boxbins[1] - boxbins[0]

        # The binedegs for halo mass binning 
        binedge = np.logspace(
        np.log(np.min(halodata["N"])),
        np.log(np.max(halodata["N"])),
        halobins + 1,
        base=math.e)

        # The location to dave the files generated
        pathsave = savelocation(params,redshift)

        # Looping for various subbox sizes
        for nw_boxsize in nw_boxsizes: 

            # Making seperate folders to save CICs
            nwboxloc = os.path.join(pathsave,str(nw_boxsize))

            # Creating the directory if it does not exist
            if not (os.path.exists(nwboxloc)):
                print("Creating directory to store data.")
                os.makedirs(nwboxloc)

            # Starting the process of cutting the box to smaller files
            for pos_start in itertools.product(boxbins,boxbins,boxbins):

                print(f"Boxsize being computed: {nw_boxsize}")
                print(f"Part of box being computed: {str(pos_start)}")
                cutside = int(dbin/nw_boxsize)
                print(f"The times the box length is cut: {cutside}")
                totalboxes = int(cutside**3) 
                print(f"The total number of boxes: {totalboxes}")


                # Filtering data for the small box
                halodata_fitr = halodata[["N","xpos","ypos","zpos"]][ (halodata["xpos"] >= pos_start[0]) & (halodata["xpos"] < (pos_start[0] + dbin)) &
                                                                 (halodata["ypos"] >= pos_start[1]) & (halodata["ypos"] < (pos_start[1] + dbin)) &
                                                                 (halodata["zpos"] >= pos_start[2]) & (halodata["zpos"] < (pos_start[2] + dbin)) ]

                partdata_fitr = partdata["xpos","ypos", "zpos"][ (partdata["xpos"] >= pos_start[0]) & (partdata["xpos"] < (pos_start[0] + dbin)) &
                                                            (partdata["ypos"] >= pos_start[1]) & (partdata["ypos"] < (pos_start[1] + dbin)) &
                                                            (partdata["zpos"] >= pos_start[2]) & (partdata["zpos"] < (pos_start[2] + dbin)) ]


                cicdata = np.zeros((totalboxes,halobins+2))

                for i in range(halobins):
                    data = halodata_fitr[["xpos","ypos","zpos"]][(binedge[i] <= halodata_fitr["N"]) & (halodata_fitr["N"] < binedge[i + 1])]
                    data = np.array([data['xpos'], data['ypos'], data['zpos']]).T
                    if len(data) == 0:
                        continue
                    cicdata[:,i] = cic(data, params, nw_boxsize, pos_start, dbin)

                partpos = np.array([partdata_fitr["xpos"], partdata_fitr["ypos"], partdata_fitr["zpos"]]).T

                cicdata[:,halobins] = cic(partpos, params, nw_boxsize, pos_start, dbin)
                cicdata[:,halobins+1] = dens_contrast(cicdata[:,halobins], params, nw_boxsize, pos_start, dbin)

                np.save(os.path.join(nwboxloc,str(pos_start)),cicdata)

            np.save(os.path.join(nwboxloc,"mass_bins"),binedge)

        totalhalos = len(halodata)
        totalpart = len(partpos)
        smhalomass = np.min(halodata["N"])  
        smhalono = np.min(halopartno)     
        lrhalomass = np.max(halodata["N"])
        lrhalono = np.max(halopartno)
        with open(os.path.join(pathsave,"mass_data.txt"),"w") as file:
            file.write('\n'.join([f"Number of halos = {totalhalos}",
                                  f"Number of particles = {totalpart}", 
                                  f"Smallest Mass of Halo = {smhalomass:.18E}", 
                                  f"No: of particles in smallest halo = {smhalono}",
                                  f"Largest Mass of halo = {lrhalomass:.18E}",
                                  f"No: of particles in largest halo =  {lrhalono}",
                                  f"Mass of a particle = {params.mass:.18E}"       ]) ) 
