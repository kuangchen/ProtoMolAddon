# USING THE NEW STRUCTURE
from MDL import *

# PHYSICAL
phys = Physical()
io = IO()
io.readPDBPos(phys, "data/alanine.pdb")
io.readPSF(phys, "data/alanine.psf")
io.readPAR(phys, "data/alanine.par")
io.readEigenvectors(phys, "data/eigVmC7eq")
phys.bc = "Vacuum"
phys.cellsize = 5
phys.exclude = "scaled1-4"
phys.temperature = 300
phys.seed = 1234


# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("badi")
ff.nonbondedForces("le")
ff2 = forces.makeForceField(phys)
ff3 = forces.makeForceField(phys)
ff2.bondedForces("badi")
ff2.nonbondedForces("le")
ff3.bondedForces("badi")
ff3.nonbondedForces("le")

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme=['NormalModeMori', 'NormalModeRelax', 'NormalModeBrownian'],
                       steps=20,
                       cyclelength=[1,1],
                       dt=4.0,
                       forcefield=[ff, ff2, ff3],
                       params={'NormalModeMori':{'firstmode':1,
                                                 'numbermodes':10,
                                                 'gamma':80,
                                                 'seed':1234,
                                                 'temperature':300,
                                                 'modeOutput':'data/test1',
                                                 'instForce':0},
                               'NormalModeRelax':{'minimlim':0.5,
                                                  'rediag':0,
                                                  'simplemin':1},
                               'NormalModeBrownian':{'firstmode':11,
                                                     'numbermodes':56,
                                                     'gamma':420,
                                                     'seed':1234,
                                                     'temperature':300,
                                                     'avForceFile':'',
                                                     'inForceFile':''}})

io.writeXYZPos(phys, 'data/NMMTest.xyz')
