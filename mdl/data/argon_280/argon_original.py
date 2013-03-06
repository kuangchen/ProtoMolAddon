# USING THE NEW STRUCTURE

from MDL import *

from propagators.examples.Leapfrog import *

ff = forcefield("charmm")
 
def physical():
    coordinates("work/untitled.pdb")
    structure("work/argon.psf")
    velocities("work/argon_eq.vel")
    parameters("work/argon.par")
    boundaryConditions("Periodic")
    temperature(106)
    #cellBasisVector1(12.5,0,0)
    #cellBasisVector1(0,12.5,0)
    #cellBasisVector1(0,0,12.5)
    #cellOrigin(0,0,0)
    cellSize(6.5)
 
def approximations():
    cutoff(ff, "LennardJones", 12.0)
    extraparameters(ff, "LennardJones", "switchon=10.0")
    algorithm(ff, "LennardJones", "Cutoff")
    switching(ff, "LennardJones", "C1")
 
def output():
    writeEnergies("a-energies.dat", 1) 
    #writeMomentum("momentum.dat", 200)

def execute():
    # RUN BBK, NUMSTEPS=20, TIMESTEP=0.5, FORCEFIELD=ff2, PLUS EXTRA PARAMETERS
    gamma = propagate("Leapfrog", 100, 10.0, ff)
	
	
