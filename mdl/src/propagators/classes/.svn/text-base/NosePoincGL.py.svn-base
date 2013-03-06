from STS import *
import Constants
import numpy

class NosePoincGL(STS):
   """
   This method is based on an extended Hamiltonian with an additional thermostat, or degree of freedom. The resulting Hamiltonian gives rise to implicit coupling between the variables, requiring implicit symplectic methods. Since this is undesirable the Generalised Leapfrog method is utilised to produce a system which can be solved explicitly. This method is described in detail in:
   
   S.D.Bond, B.Laird, and B.Leimkuhler, The Nose-Poincare method for constant temperature molecular dynamics. J.Comp.Pys.,151:114,1999.  
   """
   def init(self, phys, forces, prop):
      """
      Initialize propagator.
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
      self.bathPold = self.bathP #: Old thermostat value as system stores rescaled velocity number of Dof of system, momentum conserved so atoms*3D-3
      self.gkT = 3.0*(phys.numAtoms()-1.0)*Constants.boltzmann()*self.temp #: Product of degrees of freedom, boltzmann constant and Kelvin temperature
      self.KEtoT = 2.0 / (3.0*(phys.numAtoms()-1.0)*Constants.boltzmann()) #: Kinetic to temp. conversion
      prop.calculateForces(forces)
      self.Potnl = forces.energies.potentialEnergy(phys)  #: Potential energy
      self.h0 = self.Potnl + forces.energies.kineticEnergy(phys)  #: Initial 'internal' Hamiltonian value so that total Hamiltonian allways 0
      self.stepsdone = 0 #: Number of steps completed
      self.avTemp = 0 #: Average Kelvin temperature
      self.tempers = [] #: Holds pairs of step numbers and average temperatures
      self.Hamiltonian = [] #: Holds pairs of step numbers and total energies
 
   def half1UpdtMom(self, phys, forces, prop):
      """
      Update system velocity first 1/2 step

      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
      phys.velocities += forces.force*0.5*self.dt*phys.invmasses

   def half1UpdtbathM(self, phys, forces, prop):
      """
      Update thermal bath 'momenta' first 1/2 step

      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
      tempC = 0.5*self.dt*(self.gkT*(1.0 + numpy.log(self.bathP))-forces.energies.kineticEnergy(phys)+self.Potnl-self.h0)-self.bathM
      self.bathM = -2*tempC / (1.0 + numpy.sqrt(1.0 - tempC*self.dt/self.Q))

   def UpdtPosBathP(self, phys, prop):
      """
      Update positions and thermal bath variable full step

      @type phys: Physical
      @param phys: The physical system.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
      self.bathPold = self.bathP
      tempF = 0.5*self.dt*self.bathM/self.Q
      self.bathP *= (1.0 + tempF)/(1.0 - tempF)
      phys.positions += phys.velocities*self.bathPold*0.5*self.dt*(1.0/self.bathP+1.0/self.bathPold)

   def half2UpdtbathM(self, phys, forces,prop):
      """
      Update thermal bath 'momenta' first 1/2 step

      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
      tempS = self.bathPold/self.bathP
      tempS *= tempS
      self.bathM += 0.5*self.dt*(forces.energies.kineticEnergy(phys)*tempS-self.gkT*(1.0+numpy.log(self.bathP))-self.Potnl+self.h0-0.5*self.bathM*self.bathM/self.Q)

   def half2UpdtMom(self, phys, forces, prop):
      """
      Update system velocity second 1/2 step

      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
      tempS = self.bathPold/self.bathP
      phys.velocities *= tempS
      phys.velocities += forces.force*0.5*self.dt*phys.invmasses
                                                                                                                                                                          
   def run(self, phys, forces, prop):
      """
      Run propagator.
      1/2 step in momentum followed by 1/2 step in thermostat momentum
      followed by full step in positions, then 1/2 step in thermostat momentum
      and 1/2 step in momentum to give a time reversible method.
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """

      self.half1UpdtMom(phys, forces, prop)
      self.half1UpdtbathM(phys, forces, prop)
      self.UpdtPosBathP(phys, prop)
      prop.calculateForces(forces)
      self.Potnl = forces.energies.potentialEnergy(phys)
      self.half2UpdtbathM(phys, forces, prop)
      self.half2UpdtMom(phys, forces, prop)
      self.avTemp += forces.energies.kineticEnergy(phys)*self.KEtoT
      self.stepsdone += 1

   def finish(self, phys, forces, prop):
      """
      Finalize propagator.
      Append data to temperature and total energy Python lists.
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
      self.tempers.append([self.stepsdone,self.avTemp/self.stepsdone])
      self.Hamiltonian.append([self.stepsdone,self.bathP*(forces.energies.kineticEnergy(phys)+forces.energies.potentialEnergy(phys)+0.5*self.bathM*self.bathM/self.Q+self.gkT*numpy.log(self.bathP)-self.h0)])
      
      
name="NosePoincGL"             #: Name of propagation scheme.
parameters=('temp', 500,
            'bathM', 0.0,
            'bathP', 1.0,
            'Q', 10.0)         #: Parameters and defaults
