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
gamma = prop.propagate(scheme=['NormalModeDiagonalize', 'NormalModeLangevin', 'NormalModeMinimizer'],
                       steps=20,
                       cyclelength=[1,1],
                       dt=4.0,
                       forcefield=[ff, ff2, ff3],
                       params={'NormalModeDiagonalize':{'reDiagFrequency':100,
                                                        'minSteps':20,
                                                        'minLim':0.1,
                                                        'removeRand':1},
                               'NormalModeLangevin':{'firstmode':1,
                                                     'numbermodes':22,
                                                     'gamma':80,
                                                     'seed':1234,
                                                     'gencompnoise':0,
                                                     'temperature':300},
                               'NormalModeMinimizer':{'minimlim':0.5,
                                                      'rediag':0,
                                                      'randforce':1,
                                                      'simplemin':1,
                                                      'euFactor':0.5}})
io.writeXYZPos(phys, 'data/NMLTest.xyz')
