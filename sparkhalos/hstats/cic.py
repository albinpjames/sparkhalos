import math
import numpy as np
import os
from astropy.io import ascii

def cic_halos_old(params, redshift, nw_boxsize, totalbins):
    """This CIC function was written to calculte the cic from halos without 
    converting the astropy table to numpy arrays
    """
    cutside = params.boxsize/nw_boxsize
    totalboxes = cutside**3

    pathset = os.path.join(
        params.datadirec,
        "ProcessedData",
        "AbacusSummit_" + params.type + "_" + params.cosmo + "_" + params.intcont,
        "halos/z" + redshift,
    )

    # check the directory does not exist
    if not (os.path.exists(pathset)):
        print("Processed file does not exist")

    print("Reading file")

    # Reads the data to be used for caluclating hmf
    data = ascii.read(os.path.join(pathset, params.filename_notime + "_" + redshift + ".dat"))    
    
    x = data
    x["xpos"] = data["xpos"]/nw_boxsize
    x["ypos"] = data["ypos"]/nw_boxsize
    x["zpos"] = data["ypos"]/nw_boxsize

    # x = data[["xpos","ypos","zpos"]]/nw_boxsize
    x = x.astype(int)
    
    binedge = np.logspace(
        np.log(np.min(data["N"])),
        np.log(np.max(data["N"])),
        totalbins + 1,
        base=math.e,
    )

    dm = np.log10(binedge[1])-np.log10(binedge[0])

    boxno = x["xpos"] + cutside*x["ypos"] + cutside**2*x["zpos"]
    data.add_columns([boxno],names=["boxno"])


    #Box halo numbers
    boxdata = np.zeros((totalbins,totalboxes))

    for f in range(totalboxes):
        data_fil = data["N"][data["boxno"] == f]
        for i in range(totalbins):
            boxdata[i][f] = len(data_fil["N"][int((np.log10(data_fil["N"]) - np.log10(binedge[0]))/dm) == i])


    return boxdata

def cic_halos(data, params, redshift, nw_boxsize, totalbins):
    cutside = int(params.boxsize/nw_boxsize)
    totalboxes = int(cutside**3) 

    binedge = np.logspace(
        np.log(np.min(data["N"])),
        np.log(np.max(data["N"])),
        totalbins + 1,
        base=math.e,
    )

    dm = np.log10(binedge[1])-np.log10(binedge[0])

    data = np.array(data)
    data = np.lib.recfunctions.structured_to_unstructured(data)
    x = data[:,1:4]/nw_boxsize
    x = x.astype(int)
    
    boxes = []
    for i in range (len(x)):
        boxes.append(x[i][0] + cutside*x[i][1] + cutside**2*x[i][2])

    #Box halo numbers
    boxdata = np.zeros(totalboxes)
    for i in range(len(boxes)):
        boxdata[i] += 1
        
    return boxdata


def cic_particles(data, params, nw_boxsize, cic_method ="binned_stat", density_contrast = True):
    match cic_method:
        case "manual":
            print("Calculating CIC Using Maunal Method")
            cutside = int(params.boxsize/nw_boxsize)
            totalboxes = int(cutside**3) 
            x = data[:,0:3]/nw_boxsize
            print("converting to integer")
            x = x.astype(int)

            boxes = []
            print("caclulating boxes for particles")
            
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
                size = params.boxsize
                x_bins_dd = np.arange(0,size + nw_boxsize, nw_boxsize) 
                y_bins_dd = np.arange(0,size + nw_boxsize, nw_boxsize) 
                z_bins_dd = np.arange(0,size + nw_boxsize, nw_boxsize) 
                boxdata = binned_statistic_dd([np.array(data[:,0]),np.array(data[:,1]),np.array(data[:,2])]
                    ,values = None, statistic = 'count', 
                    bins =[x_bins_dd, y_bins_dd, z_bins_dd]).statistic

                boxdata = boxdata.ravel

    if density_contrast:
        size = params.boxsize
        # size = 200
        cell_avg = np.sum(boxdata) * (nw_boxsize**3) / (size**3)
        celldensity = (boxdata - cell_avg)/cell_avg
        oldboxdata = boxdata
        boxdata = celldensity

    return boxdata