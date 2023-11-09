"""This module defines the parameters and the functions to acess and process
the data from the abacus summit simulation. 
Data download portal: https://abacusnbody.org/ 
Documentation: https://abacussummit.readthedocs.io/en/latest/
Paper: https://academic.oup.com/mnras/article/508/3/4017/6366248?login=true
"""

import numpy as np
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
from astropy.table import Table, vstack
from astropy.io import ascii
from astropy.table import Column
import os
from pathlib import Path

from abacusnbody.data.read_abacus import read_asdf

def _readdata(clms, params, redshift):
    """This function reads the data from abacus summit simulation as a whole,
    this should be used when enough memory is available to load the complete
    data and process it together.
    """
    print("Reading the data")

    # Location of the data
    file = os.path.join(
        params.datadirec,
        "Simulations/AbacusSummit_Public_Data_Access/AbacusSummit_"
        + params.type
        + "_"
        + params.cosmo
        + "_"
        + params.intcont,
        "halos/z" + redshift,
        "halo_info/",
    )

    cat = CompaSOHaloCatalog(file, cleaned=False)  # Reads the data
    data = cat.halos[clms]  # Reads the given column and saves it to an array
    del cat  # Deleting complete loaded data to save memory

    return data  # Returning the requested colums


def _read1by1(clms, params, redshift):
    """This function reads the data from abacus summit simulation file by file,,
    this should be used when not enough memory is available to load the complete
    data and process it together.
    """
    print("Status: Reading the files")

    data = Table()
    for i in range(params.fno_s, params.fno_e + 1):
        print("Procesing file", i)
        if i < 10:
            file = os.path.join(
                params.datadirec,
                "Simulations/AbacusSummit_Public_Data_Access/AbacusSummit_"
                + params.type
                + "_"
                + params.cosmo
                + "_"
                + params.intcont,
                "halos/z" + redshift,
                "halo_info/halo_info_00" + str(i) + ".asdf",
            )
        else:
            file = os.path.join(
                params.datadirec,
                "Simulations/AbacusSummit_Public_Data_Access/AbacusSummit_"
                + params.type
                + "_"
                + params.cosmo
                + "_"
                + params.intcont,
                "halos/z" + redshift,
                "halo_info/halo_info_0" + str(i) + ".asdf",
            )

        cat = CompaSOHaloCatalog(file, cleaned=False)
        data = vstack([data, cat.halos[clms]])
        del cat

    return data


def mass_pos(params, redshift, mode="all"):
    """This function extracts the mass and position of the particles in the simulation.

    Parameters:
    params: parameters of the simulation
    redshift: Redshift being calculated
    mode: Default -all; "all" to extract all the data, "1by1" to extract data by file
    """

    print("Reading halo mass position data.")

    # Extarct mass and position
    print("Status: Extracting the mass and position of the particles")
    match mode:
        case "all":
            data = _readdata(["N", "SO_central_particle"], params, redshift)

        case "1by1":
            data = _read1by1(["N", "SO_central_particle"], params, redshift)

    # Calculating the total mass of the halo
    # print("converting to mass")
    # data["N"] = data["N"] * params.mass
    data["SO_central_particle"] += params.boxsize / 2

    data.add_columns(
        [
            data["SO_central_particle"][:, 0],
            data["SO_central_particle"][:, 1],
            data["SO_central_particle"][:, 2],
        ],
        names=["xpos", "ypos", "zpos"],
    )

    del data["SO_central_particle"]
    return data


def _readrv(params, redshift, type="field", subset="A"):
    """This function reads the field data from abacus summit simulation
    type : field or halo
    subset : A or B
    """

    # Location of the data
    path = os.path.join(
        params.datadirec,
        "Simulations/AbacusSummit_Public_Data_Access/AbacusSummit_"
        + params.type
        + "_"
        + params.cosmo
        + "_"
        + params.intcont,
        "halos/z" + redshift,
        type + "_rv_" + subset
    )

    files = Path(path).glob('*.asdf')
    for i, file in enumerate(files):
        if i == 0:
            cat = read_asdf(file, cleaned=False)
        else:
            print(f"reading {i} file")
            cat_read = read_asdf(file, cleaned=False) # Reads the data
            cat =  vstack([cat, cat_read])
            del cat_read

    cat.add_columns(
        [
            cat["pos"][:, 0],
            cat["pos"][:, 1],
            cat["pos"][:, 2],
        ],
        names=["xpos", "ypos", "zpos"],
    )
    del cat["pos"]
    del cat["vel"]

    return cat

def read_particles(params, redshift, toread):
    import itertools
    for i, pardata in enumerate(itertools.product(*toread)):
        if i == 0:
            partcl = _readrv(params, redshift,*pardata)
        else:
            partcl_read = _readrv(params, redshift,*pardata)
            partcl = vstack([partcl, partcl_read])
        print("reading data complete")

    return partcl


