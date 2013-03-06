import Vector3DBlock
import Constants
import numpy
import math

def velocityscale(phys, forces, io, steps, timestep, fg, t0):
   """
   Runs the Leapfrog method, but imposes a scaling of velocities
   after execution, to keep average kinetic energy constant.
   
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
       # Scale
       phys.velocities *= math.sqrt(t0/phys.temperature)
       # Update time
       phys.time = step*timestep       
       # Increment the step
       step = step + 1
   # Update velocities by half a step
   phys.velocities += forces.force*0.5*timestep*phys.invmasses # kick
   # Scale
   phys.velocities *= math.sqrt(t0/phys.temperature)
   # Return positions and velocities as an array of arrays.


name="velocityscale"   #: Propagator name for the factory
parameters=('T0', 300)          #: Parameters and defaults
