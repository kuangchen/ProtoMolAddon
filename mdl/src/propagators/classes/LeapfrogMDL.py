from random import *
from STS import *
from math import *
import Constants

class LeapfrogMDL(STS):
    """
    Implements the Leapfrog method.
       1. Half-timestep update of velocities.
       2. Full-timestep update of positions.
       3. Half-timestep update of velocities.
    cf. R. W. Hockney and J. W. Eastwood, Computer Simulation Using Particles.
    New York: McGraw-Hill, 1981.
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
        # Update velocities by half a step
        phys.velocities += forces.force*0.5*self.dt*phys.invmasses  # half kick
        # Update positions by a full step
        phys.positions += phys.velocities*self.dt               # drift
        # Calculate new forces with updated position/velocity
        prop.calculateForces(forces)


    def run(self, phys, forces, prop):
        """
        Run the propagator.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type forces: Forces
        @param forces: MDL Forces object
        
        @type prop: Propagator
        @param prop: MDL Propagator object
        """
        # Update velocities (full step)
        phys.velocities += forces.force*self.dt*phys.invmasses # kick
        # Update positions (full step)
        phys.positions += phys.velocities*self.dt  # drift
        # Calculate new forces
        prop.calculateForces(forces)

    def finish(self, phys, forces, prop):
        """
        Finish the propagator: do one final half-kick
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type forces: Forces
        @param forces: MDL Forces object
        
        @type prop: Propagator
        @param prop: MDL Propagator object
        """
        # Update velocities by half a step
        phys.velocities += forces.force*0.5*self.dt*phys.invmasses # kick

name="LeapfrogMDL"  #: Name of propagation scheme
parameters=() #: Parameter names and defaults
