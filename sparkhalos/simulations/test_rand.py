import numpy as np
import pandas as pd


def generate_mxyz(params, numbPoints):    

    xx = np.random.uniform(0, params.boxsize, size=(numbPoints))
    yy = np.random.uniform(0, params.boxsize, size=(numbPoints))
    zz = np.random.uniform(0, params.boxsize, size=(numbPoints))

    mass = 10**(np.random.uniform(10, 15, size=(numbPoints)))
    print(mass)
    data = np.column_stack((mass,xx,yy,zz))
    data = pd.DataFrame(data, columns = ['N','xpos','ypos','zpos'])

    return data

def generate_xyz(params, numbPoints):    

    xx = np.random.uniform(0, params.boxsize, size=(numbPoints))
    yy = np.random.uniform(0, params.boxsize, size=(numbPoints))
    zz = np.random.uniform(0, params.boxsize, size=(numbPoints))
    # zz = np.random.uniform(0, 200, size=(numbPoints))

    data = np.column_stack((xx,yy,zz))

    return data