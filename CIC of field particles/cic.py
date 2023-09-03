"""This program is used to calculate the CIC and various parameters
"""

import sys
sys.path.append('..')

from sparkhalos.simparams import SimuParams
""" Choose the simulation to be used """
from sparkhalos.simulation import test_rand as simulation
from sparkhalos.simparams.test_rand import test5000 as simparams
# from sparkhalos.simulation import abacussummit as simulation
# from sparkhalos.simparams.abacussummit import hugebase2000 as simparams

def getlocation():
    import os
    path = os.getcwd()
    from pathlib import Path
    path = Path(path).parents[1]
    return os.path.join(path, "DATA")

if __name__ == "__main__":

    """ Intilaises the simulation parameters for the simulation """
    datalocation = getlocation()
    params = SimuParams.init(datalocation, simparams)

    """ Redshifts to be computed """
    # redshifts = ["3.000","2.500","2.000"]
    redshifts = ["3.000"]

    """ The new box size & Bins for computing count in cells distribution """
    # nw_boxsizes = [25]
    nw_boxsizes = [10,20,25,50]
