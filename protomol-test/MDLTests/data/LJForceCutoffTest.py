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
phys.temperature = 310


# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.nonbondedForces("l")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                        'switching':'Cn',
                        'order':2,
                        'switchon':8.0,
                        'cutoff':10.0}

# OUTPUT
io.files = {'xyztrajforce':('data/LJForceCutoffTest.xyz',1)}

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate("LangevinImpulse", steps=0, dt=1.0, forcefield=ff, params={'temp':310, 'gamma':91})
