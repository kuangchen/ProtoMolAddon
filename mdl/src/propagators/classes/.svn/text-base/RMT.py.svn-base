#from MDL import *
from STS import *
#from umath import *
import Constants
import numpy


class RMT(STS):
   """
   This method is based on an extended Hamiltonian with 'm' additional thermostats, or degrees of freedom. The resulting Hamiltonian gives rise to implicit coupling between the variables, requiring implicit symplectic methods. Since this is undesirable a Hamiltonian splitting method is used to produce 2+m Hamiltonian systems which can be solved explicitly. This spliting is described in: http://www.nd.edu/~csweet1/campus/rmtsplit.pdf
   cf. B. J. Leimkuhlder and C. R. Sweet.  Hamiltonian Formulation for
   Recursive Multiple Thermostats in a Common Timescale.  SIAM J. Appl. Dyn.
   Sys. 4(1):187-216, 2005.
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
      prop.calculateForces(forces)
      self.Potnl = forces.energies.potentialEnergy(phys) #: Potential energy
      self.gkT = 3.0*(phys.numAtoms()-1.0)*Constants.boltzmann()*self.temp #: Number of Dof*kT, momentum conserved so number of atoms * 3D - 3 
      self.KEtoT = 2.0 / (3.0*(phys.numAtoms()-1.0)*Constants.boltzmann()) #: Convertion of kinetic energy to temperature
      self.Nf = 3.0*(phys.numAtoms()-1.0)   #: number of Dof
      self.kT = Constants.boltzmann()*self.temp    #: Boltzmann constant times Kelvin temperature
      self.h0 = self.totalEnergy(0,phys,forces)    #: Initial total energy
      self.stepsdone = 0 #: Number of steps completed
      self.avTemp = 0 #: Average Kelvin temperature
      self.tempers = [] #: Holds pairs of step numbers and average temperatures
      self.Hamiltonian = [] #: Holds pairs of step numbers and total energies

   def halfUpdtH2(self,typ,phys,forces,prop):
      """
      H2 half step.
      
      @type typ: int
      @param typ: 0 (first cycle) or 1 (second cycle)
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
#     Calc Ps
      prodS = self.prodSs(1,self.NumStats)
      ii = 1
      sumSpot = 0.0
      while (ii < self.NumStats):
         sumSpot += (self.Nf+ii)*self.kT*numpy.log(self.S[ii]) + 0.5*(1.0-self.S[ii])*(1.0-self.S[ii])/self.C[ii]
         ii += 1
      self.Ps[0] -= 0.5*self.dt*prodS*(self.Potnl+sumSpot)
      ii = 1
      while (ii < self.NumStats):
         self.Ps[ii] -= 0.5*self.dt*self.S[0]*(prodS/self.S[ii])*(self.Potnl+sumSpot+(self.Nf+ii)*self.kT-self.S[ii]*(1.0-self.S[ii])/self.C[ii])
         ii += 1
#     Calculate Momenta, test for half of cycle as using SCALED velocities
      if typ < 1:
         phys.velocities += forces.force*0.5*self.dt*phys.invmasses
#     Calc current product of s'
         self.OldProdS = prodS*self.S[0]
      else:
         tempS = self.OldProdS/(prodS*self.S[0])
         phys.velocities *= tempS
         phys.velocities += forces.force*0.5*self.dt*phys.invmasses

   def halfUpdtH3j(self,dir,prop):
      """
      H3 half step.
      
      @type dir: int
      @param dir: 0 direction j=1...M, dir=1 direction j=M...1
            
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
      kk = 1
      while (kk < self.NumStats):
         if dir > 0:
            jj = self.NumStats - kk
         else:
            jj = kk  
#     1/4 step for Ps1
         prodSlj = self.prodSs(0,jj)
         prodSgj = self.prodSs(jj+1,self.NumStats)
         a = self.dt*prodSlj/(8.0*self.Q[jj]*prodSgj)
         c = -self.Ps[jj]
         self.Ps[jj] = -2.0*c/(1.0+numpy.sqrt(1.0-4.0*a*c))
#        1/2 step for Sj
         self.Sn = self.S[jj]
         a = 0.25*self.dt*prodSlj*self.Ps[jj]/(self.Q[jj]*prodSgj)
         self.S[jj] *= (1.0+a)/(1.0-a)
#        1/2 step for Ps2-PsM
         ii = 0
         while (ii < jj):
            self.Ps[ii] -= 0.25*self.dt*((self.Sn+self.S[jj])/self.S[ii])*(0.5*prodSlj*self.Ps[jj]*self.Ps[jj]/(self.Q[jj]*prodSgj))
            ii += 1
         ii = jj+1
         while (ii < self.NumStats):
            self.Ps[ii] += 0.25*self.dt*((self.Sn+self.S[jj])/self.S[ii])*(0.5*prodSlj*self.Ps[jj]*self.Ps[jj]/(self.Q[jj]*prodSgj))
            ii += 1
#        1/4 step for Ps1
         self.Ps[jj] += 0.25*self.dt*(-0.5*prodSlj*self.Ps[jj]*self.Ps[jj]/(self.Q[jj]*prodSgj))
#
         kk += 1

####H31 half step
   def halfUpdtH31(self, prop):
      """
      H31 half step.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
#     1/4 step for Ps1
      prodS = self.prodSs(1,self.NumStats)
      a = self.dt/(8.0*self.Q[0]*prodS)
      c = -self.Ps[0]-0.25*self.dt*prodS*self.h0
      self.Ps[0] = -2.0*c/(1.0+numpy.sqrt(1.0-4.0*a*c))
#     1/2 step in S1
      self.Sn = self.S[0]
      a = 0.25*self.dt*self.Ps[0]/(self.Q[0]*prodS)
      self.S[0] *= (1.0+a)/(1.0-a)
#     1/2 step for Ps2-PsM
      ii = 1
      while (ii < self.NumStats):
         self.Ps[ii] += 0.25*self.dt*((self.Sn+self.S[0])/self.S[ii])*(0.5*self.Ps[0]*self.Ps[0]/(self.Q[0]*prodS)+prodS*self.h0)
         ii += 1
#     1/4 step for Ps1
      self.Ps[0] += 0.25*self.dt*(prodS*self.h0-0.5*self.Ps[0]*self.Ps[0]/(self.Q[0]*prodS))

####H1 full step
   def UpdtH1(self, phys, forces, prop):
      """
      H1 full step.
            
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
#     Positions
      prodS = self.prodSs(0,self.NumStats)
      tempS = self.OldProdS/prodS
      phys.positions += phys.velocities*tempS*self.dt
      #setPosition(position() + timestep(self)*velocity()*tempS)
#     Ps'
      sKinetic = forces.energies.kineticEnergy(phys)*self.OldProdS*self.OldProdS/prodS
      self.Ps[0] += (self.dt/self.S[0])*(sKinetic-prodS*self.gkT*(1.0+numpy.log(self.S[0])))
      ii = 1
      while (ii < self.NumStats):
          self.Ps[ii] += (self.dt/self.S[ii])*(sKinetic - prodS*self.gkT*numpy.log(self.S[0]))
          ii += 1


   def totalEnergy(self,typ,phys,forces):
      """
      Calculate total energy of the system.
      
      @type typ: int
      @param typ: 0 calc h0, 1 calc non time reparam Energy, 2 calc total
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.      
      """
      prodS = self.prodSs(0,self.NumStats)
      resProdS = prodS/self.S[0]
      tempH = forces.energies.kineticEnergy(phys)+self.Potnl+self.gkT*numpy.log(self.S[0])+0.5*self.Ps[0]*self.Ps[0]/(self.Q[0]*resProdS*resProdS)
      ii = 1
      while (ii < self.NumStats):
         tempH += (self.Nf+ii)*self.kT*numpy.log(self.S[ii])+0.5*(1.0-self.S[ii])*(1.0-self.S[ii])/self.C[ii]
         resProdS /= self.S[ii]
         tempH += 0.5*self.Ps[ii]*self.Ps[ii]/(self.Q[ii]*resProdS*resProdS)
         ii += 1
      if typ > 0:
         tempH -= self.h0
      if typ > 1:
         tempH *= prodS
      return(tempH)     

####Returns product of S params
   def prodSs(self,start,end):
      """
      Returns cumulative product of S parameters

      @type start: int
      @param start: Starting index

      @type end: int
      @param end: Ending index
      """
      ii = start
      prodS = 1.0
      while (ii < end):
         prodS *= self.S[ii]
         ii += 1
      return(prodS)

####Run integrator thro n steps
   def run(self, phys, forces, prop):
         """
         Run the propagator.
         Solves for half step in H2,H31,...,H3m
         then full step in H1
         followed by half steps in H3m,...,H31,H2
         so that method is time reversible.
         
         @type phys: Physical
         @param phys: The physical system.
         
         @type forces: Forces
         @param forces: MDL Forces object.
         
         @type prop: Propagator
         @param prop: MDL Propagator object.
         """

         self.halfUpdtH2(0, phys, forces, prop)
         self.halfUpdtH31(prop)
         self.halfUpdtH3j(0, prop)
         self.UpdtH1(phys, forces, prop)
         prop.calculateForces(forces)
         self.Potnl = forces.energies.potentialEnergy(phys)
         self.halfUpdtH3j(1, prop)
         self.halfUpdtH31(prop)
         self.halfUpdtH2(1, phys, forces, prop)
         #step +=  1
         self.avTemp += forces.energies.kineticEnergy(phys)*self.KEtoT
         self.stepsdone += 1

   def finish(self, phys, forces, prop):
      """
      Finalize the propagator.
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object.
      
      @type prop: Propagator
      @param prop: MDL Propagator object.
      """
#     Save av temperature
      self.tempers.append([self.stepsdone,self.avTemp/self.stepsdone])
#     Save Hamiltonian values
      self.Hamiltonian.append([self.stepsdone,self.totalEnergy(2,phys,forces)])

name = "RMT"  #: Name of the propagator.
parameters = ('temp', 500,
              'Ps', [0.0,0.0,0.0,0.0,0.0],
              'S', [1.0,1.0,1.0,1.0,1.0],
              'Q', [2.0,3.0,4.5,6.8,10.2],
              'C', [0.04,0.04,0.04,0.04,0.04],
              'NumStats', 5)  #: Parameters and defaults.

