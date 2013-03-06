import Constants

def friction(phys, forces, prop, obj):
   """
    Modify the force vector to include a frictional force.
    This can be used for incorporation of Langevin dynamics.
    
    @type phys: Physical
    @param phys: The physical system.

    @type forces: Forces
    @param forces: MDL Forces object
    
    @type prop: Propagator
    @param prop: MDL Propagator object

    @type obj: STS/MTS
    @param obj: Prototyped propagator object
   """
   targetKE = 3.0/2.0*phys.numAtoms()*Constants.boltzmann()*obj.temp
   obj.thermal = 0
   obj.bathpos = (forces.energies.kineticEnergy(phys) - targetKE)*obj.thermal/(obj.dt*phys.numAtoms())
   forces.force -= phys.velocities*obj.bathpos*phys.masses

