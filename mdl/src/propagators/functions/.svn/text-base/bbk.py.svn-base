import Constants
import Vector3DBlock
import numpy

def bbk(phys, forces, io, steps, timestep, fg, temp, gamma, seed):
   """
   Brunger-Brooks-Karplus propagation method.
   cf. A. Brunger, C. B. Brooks and M. Karplus.  Stochastic
   Boundary Conditions for Molecular Dynamics Simulations of ST2 Water.
   Chem. Phys. Lett, v. 105, pages 495-500, 1982.
   Single timestepping.
   
   @type phys: Physical
   @param phys: The physical system.

   @type forces: Forces
   @param forces: MDL Forces object.

   @type io: IO
   @param io: MDL IO object.

   @type steps: int
   @param steps: Number of steps to run.

   @type timestep: float
   @param timestep: Timestep for propagation.

   @type fg: ForceField
   @param fg: MDL force field for evaluation.

   @type temp: float
   @param temp: Kelvin temperature.

   @type gamma: float
   @param gamma: Friction constant in Langevin dynamics.

   @type seed: int
   @param seed: Random number seed
   """
   gamma = gamma*0.001/Constants.invTimeFactor()
   #fg.calculateForces(phys, forces)
   step = 0
   while (step < steps):
       forceconstant = 2*Constants.boltzmann()*temp*gamma/timestep      # assign random force
       forces.force += forces.randomForce(phys,seed)*numpy.sqrt(forceconstant*phys.masses)
       phys.velocities *= (1.0-0.5*timestep*gamma)         # first half kick
       phys.velocities += forces.force*0.5*timestep*phys.invmasses
       phys.positions += phys.velocities*timestep                # drift
       fg.calculateForces(phys, forces)
       forceconstant = 2*Constants.boltzmann()*temp*gamma/timestep      # assign random force
       forces.force += forces.randomForce(phys,seed)*numpy.sqrt(forceconstant*phys.masses)
       phys.velocities += forces.force*0.5*timestep*phys.invmasses   # second first kick
       phys.velocities *= (1.0/(1.0+0.5*timestep*gamma))
       step = step + 1


name="bbk"  #: Propagator name for the factory
parameters=("temp", 300,
            "gamma", 2,
            "seed", 1234) #: Parameter names and defaults
     
