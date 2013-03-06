# Simulation of decalanine
# Using Nose-Hoover for propagation.
from MDL import *

# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/decalanine_66/alanine.pdb")
io.readPSF(phys, "data/decalanine_66/alanine.psf")
io.readPAR(phys, "data/decalanine_66/alanine.par")
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.temperature = 300

# FORCES
forces = Forces()
ff = forces.makeForceField(phys, "charmm")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'C1',
                             'cutoff':8.0}

# OUTPUT
#io.plots = {'totalenergy':4}
io.screen = 2

# EXECUTE
prop = Propagator(phys, forces, io)
prop.propagate(scheme="NoseNVTLeapfrog", steps=2000, dt=0.5, forcefield=ff)    
