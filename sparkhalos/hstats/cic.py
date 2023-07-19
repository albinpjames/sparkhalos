from math import factorial as fact, sqrt as sq, exp, pi
from scipy.special import gamma, gammaln
import math
import numpy as np
import os
from astropy.io import ascii


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





def cic_halos(params, redshift, nw_boxsize, totalbins):
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
