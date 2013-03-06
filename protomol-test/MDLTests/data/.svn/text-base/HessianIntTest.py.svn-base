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
ff = forces.makeForceField(phys)
ff.bondedForces("badi")
ff.nonbondedForces("le")


# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate("HessianInt", steps=20, dt=0.5, forcefield=ff, params={'eigvecFile':'data/eigVmC7eq', 'eigvalFile':'data/eigEmC7eq', 'sortByAbs':1, 'textEigFile':0})
io.writeXYZPos(phys, 'data/HessianIntTest.xyz')
