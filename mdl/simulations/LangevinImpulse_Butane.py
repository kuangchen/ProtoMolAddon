# Simulation of alanine dipeptide in vacuum
# Using Langevin Impulse
# This also records execution time.

from MDL import *
from forces.HDForce import * # PYTHON FORCE CLASS
from forces.ElectrostaticForce import *
import time
start = time.time()

# PHYSICAL SYSTEM
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/alanDipeptideBlock/alanC7axial.pdb")
io.readPSF(phys, "data/alanDipeptideBlock/blockdialanine.psf")
io.readPAR(phys, "data/alanDipeptideBlock/par_all27_prot_lipid.inp")
phys.bc = "Vacuum"
phys.temperature = 300
phys.exclude = "scaled1-4"
phys.cellsize = 226.5


# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("badi")
ff.nonbondedForces("l")

# PYTHON FORCES
hd = HDForce(phys, forces, 3.14/2.0, 1, 5.0)
ff.addPythonForce(hd)

es = ElectrostaticForce(phys, forces)
ff.addPythonForce(es)

# IO
io.screen = 1
#io.plots = {'kineticenergy':2}

# EXECUTE
prop = Propagator(phys, forces, io)
prop.propagate(scheme="LangevinImpulse", steps=1000, dt=1.0, forcefield=ff)
    
stop=time.time()
print "TOTAL TIME: ", stop-start
