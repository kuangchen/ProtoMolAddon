# USING THE NEW STRUCTURE
from MDL import *

# PHYSICAL
phys = Physical()
io = IO()
io.readXYZPos(phys, "data/water_1002/water.1002.minim.xyz")
io.readPSF(phys, "data/water_1002/water.1002.psf")
io.readPAR(phys, "data/water_1002/water.1002.par")
phys.bc = "Periodic"
phys.cellsize = 6
phys.exclude = "scaled1-4"
phys.temperature = 300
phys.cB1[0] = 21.559
phys.cB2[1] = 21.559
phys.cB3[2] = 21.559
phys.cO[0] = 10.779
phys.cO[1] = 10.779
phys.cO[2] = 10.779
phys.defaultCBV = False


# FORCES
forces = Forces()
ff = forces.makeForceField(phys)
ff.bondedForces("ba")
ff.nonbondedForces("lc")
ff.params['LennardJones'] = {'algorithm':'Cutoff',
                             'switching':'Cn',
			     'cutoff':10.0,
                             'order':4,
                             'switchon':6.0}
ff.params['Coulomb'] = {'algorithm':'Cutoff',
                             'switching':'Cn',
			     'cutoff':10.0,
                             'order':4,
                             'switchon':0.0}

# OUTPUT
io.screen = 1
io.files = {'energies':('s2hmc.energies', 1)}

# EXECUTE
prop = Propagator(phys, forces, io)
gamma = prop.propagate(scheme="s2hmc", steps=500, dt=1.0, forcefield=ff,
                       params={'md_steps':100})
