import Vector3DBlock
import Constants
import numpy
# Leapfrog function accepts initial positions and velocities
# Along with the number of steps to run, the timestep and a group
# of forces
def leapfrog(phys, forces, io, steps, timestep, fg):
   """
   Implements the Leapfrog method.
      1. Half-timestep update of velocities.
      2. Full-timestep update of positions.
      3. Half-timestep update of velocities.
   cf. R. W. Hockney and J. W. Eastwood, Computer Simulation Using Particles.
   New York: McGraw-Hill, 1981.
   
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

   """
   # Update velocities by half a step
   print timestep
   phys.velocities += forces.force*0.5*timestep*phys.invmasses  # half kick
   # Update positions by a full step
   phys.positions += phys.velocities*timestep               # drift
   # Calculate new forces with updated position/velocity
   fg.calculateForces(phys, forces)
   step = 1
   # Run for the number of passed steps
   while (step < steps):
       # Run I/O
       io.run(phys, forces, step, timestep)
       # Update velocities (full step)
       phys.velocities += forces.force*timestep*phys.invmasses # kick
       # Update positions (full step)
       phys.positions += phys.velocities*timestep  # drift
       # Calculate new forces
       fg.calculateForces(phys, forces)
       # Update time
       phys.time = step*timestep       
       # Increment the step
       step = step + 1
   # Update velocities by half a step
   phys.velocities += forces.force*0.5*timestep*phys.invmasses # kick

   # Return positions and velocities as an array of arrays.
   return [phys.positions, phys.velocities]


name="leapfrog"   #: Propagator name for the factory
parameters=()     #: Parameters and defaults
