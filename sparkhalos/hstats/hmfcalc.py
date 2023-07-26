from hmf import MassFunction
from astropy.cosmology import Planck18
import os
import numpy as np
import math
from astropy.io import ascii


def hmfcalc(redshift):
    """This function calculates the hmf using the
    Tinker08 (Tinker and Kravtsov 2008) as the defulat fititng function.
    """
    print("HMF from HMF Calculator")
    redshift = float(redshift)
    hmf = MassFunction(cosmo_model=Planck18, z=redshift)
    return hmf


def dndlnm(params, redshift: float, totalbins: int):
    """Computes the hmf from the halos as per the given parameters,
    raw simulation data has to be processed before computing this.
    """
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

    # Calulates the bin edges
    print("Calculating bin edges")
    binedge = np.logspace(
        np.log(np.min(data["N"])),
        np.log(np.max(data["N"])),
        totalbins + 1,
        base=math.e,
    )

    # Calculating the number of halos in the bin
    halomf = np.zeros(totalbins)  # Variable to store calculated halo mass function
    print("Calculating number of halos in each bin")
    for i in range(totalbins):
        halomf[i] = len(
            data["N"][(binedge[i] <= data["N"]) & (data["N"] < binedge[i + 1])]
        )

    # Calculating the bin centers
    print("Calculating bin centers")
    bin_centers = []
    for i in range(len(binedge) - 1):
        bin_centers.append(np.exp(np.sqrt(np.log(binedge[i + 1]) * np.log(binedge[i]))))

    # Binsize
    binsize = []
    for i in range(len(binedge) - 1):
        binsize.append(binedge[i + 1] - binedge[i])

    # Calculating the halo mass function
    print("Calculating halo mass function")
    for i in range(totalbins):
        halomf[i] = (halomf[i] * bin_centers[i]) / (binsize[i] * params.volume)

    return halomf, bin_centers
