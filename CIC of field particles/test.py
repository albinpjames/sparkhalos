import sys
sys.path.append('..')

from sparkhalos.simparams.simparams import SimuParams
""" Choose the simulation to be used """
# from sparkhalos.simulation import test_rand as simulation
# from sparkhalos.simparams.test_rand import test5000 as simparams
from sparkhalos.simulations import abacussummit as simulation
from sparkhalos.simparams.abacussummit import hugebase2000 as simparams

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

    """ Redshifts & boxsizes to be computed """
    # redshifts = ["3.000","2.500","2.000"]
    redshifts = ["3.000"]
    # nw_boxsizes = [25]
    nw_boxsizes = [10,20,25,50]


    # cic_method = "manual"
    cic_method = "binned_stat"
    cicbins = 20
     
    """ Number of particles taken or generated """     
    particles_taken = 10000000
    take_all = True
    density_contrast = True

    for redshift in redshifts:

        """ Here the data is read to be processed """
        match params.name:
            case "test_randomnum":
                data = simulation.generate_xyz(params,particles_taken)
                print("The data is simulated for test_rand")

            case "abacussummit":
                """Read halo data"""
                halodata = simulation.mass_pos(params, redshift, mode="all")