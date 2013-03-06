from random import *
from MTS import *
from math import *
import Constants

class ImpulseMDL(MTS):
    """
    Implements the multiple-timestepping Verlet/r-RESPA method, also
    known as Impulse.  This propagator invokes an 'inner' propagator
    for a specific number of cycles per iteration, then computes its
    own forces.
    cf. H. Grubmuller, H. Heller, A. Windemuth and K. Schulten.
    Generalized Verlet Algorithm for Efficient Molecular Dyanmics
    Simulatiosn with Long-Range Interactions.  Molecular Simulation,
    vol. 6, pages 121-142, 1991.
    """
    #  init()  -------------------------------------------------------------  #
    def init(self, phys, forces, prop):
        """
        Initialize propagator: invoke the next
        propagator in the chain, and compute forces.

        @type phys: Physical
        @param phys: The physical system.

        @type forces: Forces
        @param forces: MDL Forces object.

        @type prop: Propagator
        @param prop: MDL Propagator object.
        """
        #  Use system time to seed rng.
        prop.initNext(phys, forces)
        prop.calculateForces(forces)

    #  run()  --------------------------------------------------------------  #
    def run(self, phys, forces, prop):
            """
            Run propagator.

            @type phys: Physical
            @param phys: The physical system.

            @type forces: Forces
            @param forces: MDL Forces object.

            @type prop: Propagator
            @param prop: MDL Propagator object.
            """
            phys.velocities += forces.force*0.5*self.cyclelength*prop.myPropagator.next.dt*phys.invmasses   # half kick
            # Run the next integrator in the chain, and store its results
            # in an array
            prop.runNext(phys, forces, self.cyclelength)
            prop.calculateForces(forces)
            phys.velocities += forces.force*0.5*self.cyclelength*prop.myPropagator.next.dt*phys.invmasses

    #  finish()  -----------------------------------------------------------  #
    def finish(self, phys, forces, prop ):
        """
        Finish propagator; In this case just invoke the
        finish method of the next propagator in the chain
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type forces: Forces
        @param forces: MDL Forces object.
        
        @type prop: Propagator
        @param prop: MDL Propagator object.
        """
        return
        # Invoked once at simulation finish
        #prop.finishNext(phys, forces, prop)

name="ImpulseMDL"  #: Name of propagation scheme
parameters=()  #: Tuple of parameters
