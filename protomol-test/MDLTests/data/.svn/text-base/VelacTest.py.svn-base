# USING THE NEW STRUCTURE
from MDL import *

# PHYSICAL
phys = Physical()
io = IO()
io.readXYZBinPos(phys, "data/VelacTest.coor")
io.readXYZBinVel(phys, "data/VelacTest.vel")
io.readPSF(phys, "data/VelacTest.psf")
io.readPAR(phys, "data/VelacTest.par")
phys.bc = "Periodic"
phys.cellsize = 5
phys.exclude = "1-4"
phys.temperature = 300
phys.cB1[0] = 23
phys.cB2[1] = 23
phys.cB3[2] = 23
phys.cO[0] = 6.8
phys.cO[1] = 7.3
phys.cO[2] = 8.2


# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("bad")
ff.nonbondedForces("l")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'C2',
                             'switchon':8.0,
                             'cutoff':9.0}
                             

# OUTPUT
io.files = {'dcdtrajvel':('data/VelacTest.veldcd',10)}

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate("Leapfrog", steps=30, dt=1.0, forcefield=ff)
