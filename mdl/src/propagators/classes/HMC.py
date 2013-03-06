from random import *
from MTS import *
from math import *
import Constants

class HMC(MTS):
    """
    Implements Hybrid Monte Carlo sampling.
    HMC Integrator is implemented as Multiple-Timestepping (MTS)               
    Wraps an inner integrator and takes into account
    the energy change as a result of this integration.
    The changes in position and velocity that result from this integration are
    accepted with Metropolis probability based on the resulting energy change
    P = e^{-DH/kT}, where DH is the change in energy, k is boltzmann constant
    and T is temperature. If we reject, we roll back positions and velocities.
    cf. S. Duane, A. D. Kennedy, B. J. Pendleton and D. Roweth.  Hybrid
    Monte Carlo.  Phys. Lett. B. vol. 195 pages 216-222, 1987.
    """
    #  metropolis()  -------------------------------------------------------  #
    def metropolis( self, new, curr, phys ):
        """
        Metropolis function which computes an acceptance probability
        for sets of positions depending on energy change dE.
        P = e^-dE/kt

        @type new: float
        @param new: Hypothetical new energy with flip.

        @type curr: float
        @param curr: Current energy.

        @type phys: Physical
        @param phys: The physical system.

        @rtype: int
        @return: 0 for reject, 1 for accept
        """
        deltaE = new - curr
        
        if( deltaE < 0 ):
            return 1

        acceptProb = exp( -deltaE / ( Constants.boltzmann() * phys.getTemperature() ) )

        randNum = random()

        if( randNum < acceptProb ):
            print "\n****    Move accepted\n"
            return 1
        else:
            print "\n****    Move rejected\n"
            return 0


    #  init()  -------------------------------------------------------------  #
    def init(self, phys, forces, prop):
        """
        Initialize propagator: seed the generator, invoke the next
        propagator in the chain, and compute forces.

        @type phys: Physical
        @param phys: The physical system.

        @type forces: Forces
        @param forces: MDL Forces object.

        @type prop: Propagator
        @param prop: MDL Propagator object.
        """
        #  Use system time to seed rng.
        seed() 
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

            #  Save current positions and velocities in case we reject.
            currPos = phys.positions.copy()
            currVel = phys.velocities.copy()

            currEnergy = forces.energies.totalEnergy(phys)

            prop.runNext( phys, forces, self.cyclelength )

            newEnergy = forces.energies.totalEnergy(phys)

            accept = self.metropolis( newEnergy, currEnergy, phys )

            if( accept ):
                currPos = phys.positions.copy()
                currVel = phys.velocities.copy()

            else:
                for ii in range(0, phys.numAtoms()*3):
                    phys.positions[ii] = currPos[ii]
                    phys.velocities[ii] = currVel[ii]


    #  finish()  -----------------------------------------------------------  #
    def finish(self, phys, forces, prop ):
        """
        Finish propagator; in this case just invoke the
        finish method of the next propagator in the chain
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type forces: Forces
        @param forces: MDL Forces object.
        
        @type prop: Propagator
        @param prop: MDL Propagator object.
        """
        # Invoked once at simulation finish
        prop.finishNext(phys, forces, prop)

name="HMC"  #: Name of propagation scheme
parameters=()  #: Tuple of parameters
