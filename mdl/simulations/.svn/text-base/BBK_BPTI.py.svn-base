#from Physical import *
#from Forces import *
#from Propagator import *
#from IO import *
#from ForceField import *
from MDL import *

import time
start = time.time()

# PHYSICAL SYSTEM
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/bpti_water_1101/bpti.pdb")
io.readPSF(phys, "data/bpti_water_1101/bpti.psf")
io.readPAR(phys, "data/bpti_water_1101/bpti.par")
phys.bc = "Periodic"
phys.exclude = "scaled1-4"
phys.cellsize = 6.5
phys.temperature = 300

# FORCES
forces = Forces()
ff = forces.makeForceField(phys, "charmm")
#ff.params['Coulomb'] = {'algorithm':'PME',
#                        'switching':'Cutoff',
#                        'gridsize':15,
#                        'cutoff':12}

ff.params['Coulomb'] = {'algorithm':'MultiGrid',
                        'switching':'Cutoff',
                        'gridsize':15,
                        'smoothdist':1.0,
                        'levels':3}

# OUTPUT
io.screen = 1
#io.files = {'gui':('MDL',1)}

# PROPAGATION
prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme="BBKMDL", steps=10000, dt=1.0, forcefield=ff)
