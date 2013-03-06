# Simulation of alanine dipeptide in vacuum
# Applies normal mode restraints.
from MDL import *

phys = Physical()
io = IO()
io.readPDBPos(phys, "data/alanDipeptideVac/minC7eq.pdb")
io.readPSF(phys, "data/alanDipeptideVac/alan_mineq.psf")
io.readPAR(phys, "data/alanDipeptideVac/par_all27_prot_lipid.inp")
io.readEigenvectors(phys, "data/alanDipeptideVac/eigVmC7eq")
phys.bc = "Vacuum"
phys.temperature = 300
phys.cellsize = 5.0
phys.exclude = "scaled1-4"
phys.seed = 1234



forces = Forces()
ff = forces.makeForceField(phys, "charmm")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'C2',
                             'cutoff':12,
                             'switchon':9}
ff.params['Coulomb'] = {'algorithm':'SCPISM',
                        'switching':'Cutoff',
                        'cutoff':12}


ff2 = forces.makeForceField(phys, "charmm")
ff2.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'C2',
                             'cutoff':12,
                             'switchon':9}
ff2.params['Coulomb'] = {'algorithm':'SCPISM',
                        'switching':'Cutoff',
                        'cutoff':12}


ff3 = forces.makeForceField(phys, "charmm")
ff3.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'C2',
                             'cutoff':12,
                             'switchon':9}
ff3.params['Coulomb'] = {'algorithm':'SCPISM',
                        'switching':'Cutoff',
                        'cutoff':12}

io.screen = 1
io.files = {'gui':('MDL', 1)}

prop = Propagator(phys, forces, io)

prop.propagate(scheme=['NormalModeDiagonalize', 'NormalModeLangevin', 'NormalModeMinimizer'],
                       steps=2,
                       cyclelength=[1,1],
                       dt=4.0,
                       forcefield=[ff, ff2, ff3],
                       params={'NormalModeDiagonalize':{'reDiagFrequency':100,
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


