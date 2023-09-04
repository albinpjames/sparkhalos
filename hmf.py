from sparkhalos.simulprocess import abacussummit
from sparkhalos.simulprocess.simparams import SimuParams
from sparkhalos.hstats.hmfcalc import dndlnm, hmfcalc
from sparkhalos.hstats.cic import cic_halos


datalocation = "/mnt/dark/Projects/3.CUSAT/Data"
params = SimuParams.init(datalocation, abacussummit.simparams, type= "small", intcont= "ph3000", boxsize=500)

# Redshifts to be computed
# redshifts = ["3.000","2.000","1.100","0.500","0.200"]
redshifts = ["3.000"]
totalbins = 50

for redshift in redshifts:
    hmf = hmfcalc(redshift)
    data = mass_pos(params, redshift, mode="all")
    halomf, bin_centers = dndlnm(data, params, redshift, totalbins)

# Ploting HMF
import matplotlib.pyplot as plt
plt.scatter(bin_centers, halomf, color = '#2ab0ff' , s=10)
plt.loglog(hmf.m,hmf.dndlnm, linewidth=0.5, color ='k')