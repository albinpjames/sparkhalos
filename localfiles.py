import os

def getlocation():
    path = os.getcwd()
    from pathlib import Path
    path = Path(path).parents[1]
    return os.path.join(path, "DATA")

def savelocation(params,redshift):
    return os.path.join(
            params.datadirec,
            "ProcessedData",
            params.name,
            params.type + "_" + str(params.boxsize) + "HMpc_" + params.cosmo + "_" + params.intcont,
            "halos",
            "z" + redshift,
            "cic")
