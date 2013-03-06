# USING THE NEW STRUCTURE
from MDL import *

# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/PDBTest.in.pdb")
io.readPSF(phys, "data/PDBTest.psf")
io.readPAR(phys, "data/PDBTest.par")
phys.bc = "Vacuum"
phys.cellsize = 5
phys.exclude = "scaled1-4"
phys.temperature = 300
phys.seed = 1234


# FORCES
forces = Forces()
ff = forces.makeForceField(phys, "charmm")

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate("Leapfrog", steps=20, dt=0.5, forcefield=ff)
io.writePDBPos(phys, 'data/PDBTest.pdb')
