import numpy as np
from astropy.table import Table


def generate_mxyz(params, numbPoints):    

    xx = np.random.uniform(0, params.boxsize, size=(numbPoints))
    yy = np.random.uniform(0, params.boxsize, size=(numbPoints))
    zz = np.random.uniform(0, params.boxsize, size=(numbPoints))

    mass = 10**(np.random.uniform(10, 15, size=(numbPoints)))/params.mass
    n = mass.astype(int)
    data = np.column_stack((n,xx,yy,zz))
    data = Table(data, names= ['N','xpos','ypos','zpos'])

    return data

def generate_xyz(params, numbPoints):    

    xx = np.random.uniform(0, params.boxsize, size=(numbPoints))
    yy = np.random.uniform(0, params.boxsize, size=(numbPoints))
    zz = np.random.uniform(0, params.boxsize, size=(numbPoints))
    # zz = np.random.uniform(0, 200, size=(numbPoints))

    data = np.column_stack((xx,yy,zz))
    data = Table(data, names= ['xpos','ypos','zpos'])

    return data