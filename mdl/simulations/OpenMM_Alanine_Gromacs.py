# USING THE NEW STRUCTURE
from MDL import *

# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/alanine_gromacs/alanine.pdb")
io.readGromacs(phys, "data/alanine_gromacs/alanine.top", "data/alanine_gromacs/ffamber96/")
phys.bc = "Vacuum"
phys.temperature = 10


# FORCES
forces = Forces()
ff = forces.makeForceField(phys)

# OUTPUT
io.screen = 2

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate("OpenMM", steps=2000, dt=0.5, forcefield=ff)
