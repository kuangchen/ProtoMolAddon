# Simulation of united-atom butane
# Uses Nose-Poincare (Python prototyped)
# to propagate the system
from MDL import *

# PHYSICAL SYSTEM
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/UA_butane/UA_butane.pdb")
io.readPSF(phys, "data/UA_butane/UA_butane.psf")
io.readPAR(phys, "data/UA_butane/UA_butane.par")
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.temperature = 300

# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("ba")

# OUTPUT
#io.plots = {'totalenergy':2000}
io.screen = 4000

# PROPAGATION
prop = Propagator(phys, forces, io)
prop.propagate(scheme="NosePoincGL", steps=8000, dt=0.1, forcefield=ff)

    
