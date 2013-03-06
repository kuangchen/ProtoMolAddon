from random import *
from STS import *
from math import *
import Constants
import _TopologyUtilities
import numpy

class BBKMDL(STS):
    """
    Implements BBK propagation.
    cf. A. Brunger, C. B. Brooks and M. Karplus.  Stochastic
    Boundary Conditions for Molecular Dynamics Simulations of ST2 Water.
    Chem. Phys. Lett, v. 105, pages 495-500, 1982.
    """
    def init(self, phys, forces, prop):
       """
       Set and initialize the propagator.
       
       @type phys: Physical
       @param phys: The physical system.
       
       @type forces: Forces
       @param forces: MDL Forces object

       @type prop: Propagator
       @param prop: MDL Propagator object
       """
       self.gamma = self.gamma*0.001/Constants.invTimeFactor()  #: Dissipative factor
       prop.calculateForces(forces)

    def run(self, phys, forces, prop):
       """
       Execute the propagator.
       
       @type phys: Physical
       @param phys: The physical system.
       
       @type forces: Forces
       @param forces: MDL Forces object

       @type prop: Propagator
       @param prop: MDL Propagator object
       """
       forceconstant = 2*Constants.boltzmann()*self.temp*self.gamma/self.dt      # assign random force
       forces.force += forces.randomForce(phys,self.seed)*numpy.sqrt(forceconstant*phys.masses)
       phys.velocities *= (1.0-0.5*self.dt*self.gamma)         # first half kick
       phys.velocities += forces.force*0.5*self.dt*phys.invmasses
       phys.positions += phys.velocities*self.dt                # drift
       prop.calculateForces(forces)
       forceconstant = 2*Constants.boltzmann()*self.temp*self.gamma/self.dt      # assign random force
       forces.force += forces.randomForce(phys,self.seed)*numpy.sqrt(forceconstant*phys.masses)
       phys.velocities += forces.force*0.5*self.dt*phys.invmasses   # second first kick
       phys.velocities *= (1.0/(1.0+0.5*self.dt*self.gamma))

name="BBKMDL"  #: Name of propagation scheme
parameters=("temp", 300,
            "gamma", 2,
            "seed", 1234) #: Parameter names and defaults
