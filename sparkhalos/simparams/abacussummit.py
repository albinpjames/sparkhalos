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

hugebase2000 = {
    "name": "abacussummit",  # The name of the simulation
    "type": "hugebase",  # The type of simulation (base, small for abacussumit)
    "cosmo": "c000",  # The cosmology used in the simulation
    "intcont": "ph000",  # The initial condition used for the simulation
    "boxsize": 2000,  # The box size of the simulation
    "fno_s": 0,  # If data fragamneted to multiple files, starting file number
    "fno_e": 33,  # If data fragamneted to multiple files, ending file number
    "mass": 2.109081520453063 * 10**9,  # The mass of the particles in the simulation
}

hugebase500 = {
    "name": "abacussummit",  # The name of the simulation
    "type": "hugebase",  # The type of simulation (base, small for abacussumit)
    "cosmo": "c000",  # The cosmology used in the simulation
    "intcont": "ph000",  # The initial condition used for the simulation
    "boxsize": 500,  # The box size of the simulation
    "fno_s": 0,  # If data fragamneted to multiple files, starting file number
    "fno_e": 33,  # If data fragamneted to multiple files, ending file number
    "mass": 2.109081520453063 * 10**9,  # The mass of the particles in the simulation
}