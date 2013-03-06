# SIMULATION OF 280-ARGON MOLECULES
# USING SELF-CONSISTENT LEAPFROG
from MDL import *

# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/argon_280/untitled.pdb")
io.readPSF(phys, "data/argon_280/argon.psf")
io.readPAR(phys, "data/argon_280/argon.par")
phys.bc = "Periodic"
phys.cellsize = 6.5
phys.temperature = 300

# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.nonbondedForces("l")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                                    'switching':'C1',
                                    'cutoff':8.0}

# OUTPUT
#io.plots = {'totalenergy':4}
io.screen = 2

# EXECUTE
prop = Propagator(phys, forces, io)
prop.propagate(scheme="DMDLeapfrog", steps=200, dt=0.5, forcefield=ff)

