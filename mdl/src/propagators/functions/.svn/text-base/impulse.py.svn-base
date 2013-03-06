import Constants
import Vector3DBlock

def impulse(phys, forces, io, steps, cyclelength, fg, nextinteg, *args):
   """
   Verlet/r-RESPA propagation method.
   Implements the multiple-timestepping Verlet/r-RESPA method, also
   known as Impulse.  This propagator invokes an 'inner' propagator
   for a specific number of cycles per iteration, then computes its
   own forces.
   cf. H. Grubmuller, H. Heller, A. Windemuth and K. Schulten.
   Generalized Verlet Algorithm for Efficient Molecular Dyanmics
   Simulatiosn with Long-Range Interactions.  Molecular Simulation,
   vol. 6, pages 121-142, 1991.
   
   @type phys: Physical
   @param phys: The physical system.

   @type forces: Forces
   @param forces: MDL Forces object.

   @type io: IO
   @param io: MDL IO object.

   @type steps: int
   @param steps: Number of steps to run.

   @type cyclelength: float
   @param cyclelength: Number of iterations of inner method.

   @type fg: ForceField
   @param fg: MDL force field for evaluation.

   @type nextinteg: function handle
   @param nextinteg: Method handle for next propagator in the chain

   @type args: tuple
   @param args: Parameters for the next propagator in the chain

   """
   # Calculate initial forces
   step = 0
   # For all steps
   timestep = cyclelength*args[0]
   args2 = (args[0]*Constants.invTimeFactor(),)+args[1:len(args)]

   while (step < steps):
      # Update velocities by half a step
      phys.velocities += forces.force*0.5*timestep*phys.invmasses   # half kick
      # Run the next integrator in the chain, and store its results
      # in an array
      nextinteg(phys, forces, io, cyclelength/Constants.invTimeFactor(), *args2)
      # Calculate new forces
      fg.calculateForces(phys, forces)
      # Update velocities by another half step
      phys.velocities += forces.force*0.5*timestep*phys.invmasses   # half kick
      step = step + 1

name="impulse"  #: Propagator name for the factory
parameters=()   #: Parameter names and defaults
     
