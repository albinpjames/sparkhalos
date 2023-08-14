simparams = {
    "name": "test_randomnum",  # The name of the simulation
    "type": "test_rand",  # The type of simulation (base, small for abacussumit)
    "cosmo": "test_rand",  # The cosmology used in the simulation
    "intcont": "test_rand",  # The initial condition used for the simulation
    "boxsize": 500,  # The box size of the simulation
    "fno_s": None,  # If data fragamneted to multiple files, starting file number
    "fno_e": None,  # If data fragamneted to multiple files, ending file number
    "mass": 2.109081520453063 * 10**9,  # The mass of the particles in the simulation
}
"""Here the parameters of the abacus summit simulation are defined.

name (str): 
    Default - abacussumit; 
    The name of the simulation
type (str):  
    Default - base; 
    The type of simulation (base, small for abacussumit)
cosmo (str): 
    Default - c000; 
    The cosmology used in the simulation
intcont (str): 
    Default - ph000; 
    The initial condition used for the simulation
boxsize (int): 
    Default - 2000;  
    The box size of the simulation
fno_s (int): 
    Default - 0; 
    The starting file number, if data fragamneted to multiple files.
fno_e (int): 
    Default - 33; 
    The ending file number, if data fragamneted to multiple files.
mass (float): 
    Default - 2.109081520453063 * 10**9; 
    The mass of the particles in the simulation
"""

import numpy as np


def generate_mxyz(params, numbPoints):    

    xx = np.random.uniform(0, params.boxsize, size=(numbPoints))
    yy = np.random.uniform(0, params.boxsize, size=(numbPoints))
    zz = np.random.uniform(0, params.boxsize, size=(numbPoints))

    mass = 10**(np.random.uniform(10, 15, size=(numbPoints)))
    print(mass)
    data = np.column_stack((mass,xx,yy,zz))

    return data

def generate_xyz(params, numbPoints):    

    xx = np.random.uniform(0, params.boxsize, size=(numbPoints))
    yy = np.random.uniform(0, params.boxsize, size=(numbPoints))
    zz = np.random.uniform(0, params.boxsize, size=(numbPoints))
    # zz = np.random.uniform(0, 200, size=(numbPoints))

    data = np.column_stack((xx,yy,zz))

    return data