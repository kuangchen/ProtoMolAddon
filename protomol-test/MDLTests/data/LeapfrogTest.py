# USING THE NEW STRUCTURE
from MDL import *

# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/alanine.pdb")
io.readPSF(phys, "data/alanine.psf")
io.readPAR(phys, "data/alanine.par")
phys.bc = "Vacuum"
phys.cellsize = 5
phys.exclude = "scaled1-4"
phys.temperature = 300


# FORCES
forces = Forces()
ff = forces.makeForceField(phys, "charmm")

# EXECUTE
prop = Propagator(phys, forces, io)
prop.propagate("Leapfrog", steps=20, dt=0.5, forcefield=ff)
io.writeXYZPos(phys, 'data/LeapfrogTest.xyz')
