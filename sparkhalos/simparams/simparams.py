"""This module intaliases the variables and various parameters to be used 
while processing the simulation data. Here the variables (abacus summit has
been used as the template) are generalised for the simulations the values are
defined in the various simulation processing modules as dictonaries and are 
intialised when called.
"""

from datetime import datetime

now = datetime.now()
timestamp = now.strftime("%d-%m-%Y_%H-%M-%S")
# timestamp = datetime.isoformat(now)

from dataclasses import dataclass
import os


@dataclass
class SimuParams:
    """This class represents the parameters and directory structure of various
    simulations to be processed.
    """

    # Basic parameters of the simulation
    name: str  # The name of the simulation
    type: str  # The type of simulation (base, small for abacussumit)
    cosmo: str  # The cosmology used in the simulation
    intcont: str  # The initial condition used for the simulation
    boxsize: int  # The box size of the simulation
    # fno_s: int  # If data fragamneted to multiple files, starting file number
    # fno_e: int  # If data fragamneted to multiple files, ending file number
    mass: float  # The mass of the particles in the simulation

    volume: int # The total volume of the simulation
    # Directory structure and naming
    datadirec: str  # Directory where the simulation data is stored
    filename: str  # The name used to save the processed data
    filename_notime: str  # The name used to save the processed data without time

    @classmethod
    def init(
        cls,
        datadirectory : str,
        data,
        name : str = None,
        type : str = None,
        cosmo : str = None,
        intcont : str = None,
        boxsize : int = None,
        fno_s : int = None,
        fno_e : int = None,
        mass : float = None,
    ):
        return cls(
            name = name if name is not None else data["name"],
            type = type if type is not None else data["type"],
            cosmo = cosmo if cosmo is not None else data["cosmo"],
            intcont = intcont if intcont is not None else data["intcont"],
            boxsize = boxsize if boxsize is not None else data["boxsize"],
            # fno_s = fno_s if fno_s is not None else data["fno_s"],
            # fno_e = fno_e if fno_e is not None else data["fno_e"],
            mass = mass if mass is not None else data["mass"],
            
            volume = boxsize**3 if boxsize is not None else data["boxsize"]**3,
            datadirec = datadirectory,
            filename = str(data["boxsize"]) + "HMpc_" + data["cosmo"] + '_' + data["intcont"] + '_' + str(timestamp),
            filename_notime = str(data["boxsize"]) + "HMpc_" + data["cosmo"] + '_' + data["intcont"] 
        )
