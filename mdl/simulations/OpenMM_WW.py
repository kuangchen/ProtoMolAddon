# Recursive Multiple Thermostat technique
from MDL import *

import time
start = time.time()
# PHYSICAL SYSTEM
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/wwd/topol.pdb")
io.readGromacs(phys, "data/wwd/topol.top", "data/wwd", "data/params.agb")
#io.readPSF(phys, "data/calmodulin/calmodulin.psf")
#io.readPAR(phys, "data/calmodulin/calmodulin.par")
phys.bc = "Vacuum"
phys.cellsize = 6.5
phys.temperature = 10
phys.exclude = "scaled1-4"

# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
#ff = forces.makeForceField(phys, "charmm")
ff.gbsa = True

# OUTPUT
io.screen = 1000

# PROPAGATION
prop = Propagator(phys, forces, io)
prop.propagate(scheme="OpenMM", steps=10000, dt=0.5, forcefield=ff,
               params={'HarmonicBondForce':True,
	               'HarmonicAngleForce':True,
		       'RBDihedralForce':True,
		       'PeriodicTorsion':True,
		       'NonbondedForce':True})

stop = time.time()
print stop-start

