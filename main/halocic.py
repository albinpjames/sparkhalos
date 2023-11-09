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
    nw_boxsize = 25
    boxdata = cic_halos(params, redshift, nw_boxsize, totalbins)

import matplotlib.pyplot as plt
bins=10
plt.hist(np.sum(boxdata, axis=0), bins=bins, density=True, facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha = 0.15)