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
ff.nonbondedForces("c")
ff.params['Coulomb'] = {'algorithm':'SCPISM',
                        'BornOnly':True,
                        'switching':'Cutoff',
                        'cutoff':5,
                        'bornswitch':3}


# OUTPUT
io.files = {'xyztrajforce':('data/BornForceTest.xyz',1)}

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate("Leapfrog", steps=0, dt=0.5, forcefield=ff)
