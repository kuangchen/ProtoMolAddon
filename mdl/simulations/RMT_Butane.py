# Recursive Multiple Thermostat technique
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
ff = forces.makeForceField(phys, "charmm")

# OUTPUT
io.screen = 2

propFactory.query("RMT")
# PROPAGATION
prop = Propagator(phys, forces, io)
prop.propagate(scheme="RMT", steps=200, dt=0.5, forcefield=ff)



