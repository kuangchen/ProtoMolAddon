import MathUtilities
from PropagatorFactory import *
import time
import Constants
import TopologyUtilities
import numpy
import ProtoMolApp
import copy
import ModifierShake

class Propagator:
   """
   Provides functionality to propagate a system with time.
   """
   def __init__(self, phys, forces, io):
      #####################################################################################
      # USER-ACCESSIBLE STRUCTURES
      self.myPropagator = 0                  #: Propagator object
      self.myStep = 0                        #: Current simulation step
      self.myTimestep = 0                    #: Propagator timestep
      self.myLevel = 0                       #: Number of levels
      #####################################################################################

      phys.build()

      self.phys = phys #: Physical object 
      self.forces = forces #: Forces object
      self.io = io  #: IO object
      io.phys = phys

      
   def reset(self):
      """
      Reset the state of the Propagator object.
      """
      self.myPropagator = 0
      self.myStep = 0
      self.myTimestep = 0
      self.myLevel = 0

   def isMDL(self, integ):
      """
      Determine whether or not a propagator has been coded in MDL.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @rtype: boolean
      @return: True if the passed propagator is prototyped in pure Python.
      """
      return (hasattr(integ, "prerunmodifiers"))


   def addPreInitModifier(self, integ, modifier):
      """
      Add a modifier to execute before propagator initialization.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.preinitmodifiers.append(modifier)


   def addPostInitModifier(self, integ, modifier):
      """
      Add a modifier to execute after propagator initialization.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.postinitmodifiers.append(modifier)


   def addPreRunModifier(self, integ, modifier):
      """
      Add a modifier to execute before propagator execution.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.prerunmodifiers.append(modifier)


   def addPostRunModifier(self, integ, modifier):
      """
      Add a modifier to execute after propagator execution.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.postrunmodifiers.append(modifier)



   def addPreForceModifier(self, integ, modifier):
      """
      Add a modifier to execute before force calculation.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.preforcemodifiers.append(modifier)


   def addPostForceModifier(self, integ, modifier):
      """
      Add a modifier to execute after force calculation.

      @type integ: STS/MTS
      @param integ: MDL propagator object (STS or MTS)

      @type modifier: function
      @param modifier: Routine which alters propagator behavior.
      """
      integ.postforcemodifiers.append(modifier)


   # RUN A LIST OF PASSED MODIFIERS ON THE PASSED PROPAGATOR
   def runModifiers(self, modifiers, phys, forces, prop, integ):
      """
      Run modifiers of a propagator

      @type modifiers: list of functions
      @param modifiers: A set of routines which alternates propagator behavior

      @type phys: Physical
      @param phys: MDL Physical object

      @type forces: Forces
      @param forces: MDL Forces object

      @type prop: Propagator
      @param prop: MDL Propagator object
      
      @type integ: object
      @param integ: MDL propagator object (STS or MTS)
      """
      #integ.postf
      for ii in range(0, modifiers.__len__()):
         modifiers[ii](phys, forces, prop, integ)

   def timestep(self, integ):
      """
      Return the timestep of a propagator, scaled accordingly

      @type integ: object
      @param integ: MDL propagator object (STS or MTS)

      @rtype: float
      @return: The timestep (dt) of a propagator
      """
      return integ.getTimestep() * Constants.invTimeFactor()

   def calculateForces(self, forces):
        """
        Calculate forces and update the atomic force vector.
      
        @type forces: Forces
        @param forces: MDL Forces object
        """
        for ii in range(0, self.myPropagator.preforcemodifiers.__len__()):
           self.myPropagator.preforcemodifiers[ii](self.phys, forces, self, self.myPropagator)
        self.myPropagator.calculateForces()
        forces.forcevec = self.myPropagator.getForces()
        for ii in range(0, self.myPropagator.postforcemodifiers.__len__()):
           self.myPropagator.postforcemodifiers[ii](self.phys, forces, self, self.myPropagator)

   def initNext(self, phys, forces):
      """
      For multiple timestepping, initialize the next propagator
      in the chain.

      @type phys: Physical
      @param phys: MDL Physical object
      
      @type forces: Forces
      @param forces: MDL Forces object

      """
      tempI = self.myPropagator
      self.myLevel += 1
      setPropagator(self, phys, forces, self.myPropagator.next, levelswitch=True)
      self.myPropagator = tempI
      self.myLevel -= 1

   def runNext(self, phys, forces, cL):
      """
      For multiple timestepping, execute the next propagator
      in the chain.

      @type phys: Physical
      @param phys: MDL Physical object
      
      @type forces: Forces
      @param forces: MDL Forces object

      @type cL: int
      @param cL: Cycle length (number of times to execute the
                 inner propagator)
      """
      tempI = self.myPropagator
      self.myPropagator = self.myPropagator.next
      self.myLevel = self.myLevel + 1
      executePropagator(self, phys, forces, self.io, cL)
      self.myLevel = self.myLevel - 1
      self.myPropagator = tempI


   def finishNext(self, phys, forces, prop):
      """
      For multiple timestepping, finish the next propagator
      in the chain.

      @type phys: Physical
      @param phys: MDL Physical object
      
      @type forces: Forces
      @param forces: MDL Forces object

      @type prop: Propagator
      @param prop: MDL Propagator object

      """
      tempI = self.myPropagator
      self.myPropagator = self.myPropagator.next
      if (self.isMDL(self.myPropagator)):
         self.myPropagator.finish(phys, forces, prop)
      self.myPropagator = tempI

   def deepCopyForce(self, ff):
      p = ff.params
      if (ff.charmm):
         ff2 = self.forces.makeForceField(self.phys, "charmm")
      else:
         ff2 = self.forces.makeForceField(self.phys)
      for ii in ff.forcetypes:
         ff2.forcetypes.append(ii)
      ff2.params = p
      ff2.build()
      #self.forces.removeForceField(ff)
      import ForceGroup
      ForceGroup._swig_setattr_nondynamic(ff2, ForceGroup.ForceGroup, "thisown", 0)
      return ff2

   # PROPAGATE THE SYSTEM
   # USE METHOD "name"
   # arg1 = NUMBER OF STEPS
   # arg2 = TIMESTEP
   # arg3 = ForceField STRUCTURE
   # args = OPTIONAL EXTRA ARGUMENTS AS TUPLES
   def propagate(self, scheme="Leapfrog", steps=0, cyclelength=-1, dt=0.1, forcefield=[], params={}):
       """
       Propagate the system.

       @type scheme: string
       @param scheme: Name of the propagator to use.
       
       @type steps: int
       @param steps: Number of steps for execution.

       @type cyclelength: int
       @param cyclelength: Cyclelength for MTS propagation (-1 is STS).
       
       @type dt: float
       @param dt: Timestep.

       @type forcefield: ForceField
       @param forcefield: MDL ForceField object.

       @type params: dictionary
       @param params: Extra parameters unique for this propagation scheme.
                     (This could be empty).

       """
       self.myTimestep = dt
       chain = ()
       if (cyclelength != -1):  # MTS
          if (str(type(cyclelength))[7:11] == 'list'): # LIST, MANY LEVELS
             levels = len(cyclelength) + 1
             outertime = cyclelength[0]
          else: # ONE CYCLELENGTH = 2 LEVELS OF PROPAGATION
             levels = 2
             outertime = cyclelength

          if (str(type(scheme))[7:11] == 'list'): # LIST, MANY LEVELS
             outerscheme = scheme[0]
          else: # ONE CYCLELENGTH = 2 LEVELS OF PROPAGATION
             outerscheme = scheme


          # THE NUMBER OF FORCEFIELDS PROVIDED MUST EQUAL THE NUMBER
          # OF PROPAGATION LEVELS
          if (len(forcefield) != levels):
             print "[MDL] Error in propagate(): ", levels, " levels of propagation with ", len(forcefield), " force fields."
          outerforcefield = self.deepCopyForce(forcefield[0])

          if (str(type(scheme))[7:11] != 'list'):
             chain += (params,)
          else:
             if (params.has_key(outerscheme)):
                 chain += (params[outerscheme],)
             else:
                 chain += ({},)
	  for i in range(1, levels):
             if (str(type(scheme))[7:11] == 'list' and i < len(scheme)):
                chain += (scheme[i],)
             if (str(type(cyclelength))[7:11] == 'list' and i < len(cyclelength)):
                chain += (cyclelength[i],)
             else:
                chain += (dt,)
	     chain += (self.deepCopyForce(forcefield[i]),)
             if params.has_key(scheme[i]):
                chain += (params[scheme[i]],)
             else:
                chain += ({},)
       else: #STS
          outertime = dt
          outerscheme = scheme
          outerforcefield = self.deepCopyForce(forcefield)
          chain += (params,)
       # Build force fields.
       # Tricky, because we could be dealing with
       # a single object or list.
       if ((str(type(forcefield)))[7:len(str(type(forcefield)))-2] == 'list'):
          for ff in forcefield:
             if (ff.dirty): ff.build()
	     if (ff.gbsa): 
	        self.phys.myTop.implicitSolvent = 2
		self.phys.myTop.doGBSAOpenMM = 1
		self.phys.build()
       else:
          if (forcefield.dirty): 
              forcefield.build()
	  if (forcefield.gbsa):
	      self.phys.myTop.implicitSolvent = 2
              self.phys.myTop.doGBSAOpenMM = 1
	      self.phys.build()
       if (self.io.dirty):
          self.io.build()
        

       if (propFactory.getType(outerscheme) == "method"):
          # Calculate the forces, store them in force.
          if (not hasattr(self.phys, 'app')):
             self.phys.app = ProtoMolApp.ProtoMolApp()
             self.phys.app.makeApp(self.phys.myTop, self.phys.posvec, self.phys.velvec, self.forces.energies, dt)
          self.phys.updateCOM_Momenta()
          outerforcefield.calculateForces(self.phys, self.forces)
          self.io.run(self.phys, self.forces, 0, outertime)
          self.io.myProp = self
          for ii in range(1, steps+1):
             self.phys.app.energies.clear()
             self.forces.energies.initialize(self.phys)
             propFactory.create(1, outerscheme, self.phys, self.forces, self.io, 1, outertime*Constants.invTimeFactor(), outerforcefield, *chain)
             self.phys.time = ii*outertime*Constants.invTimeFactor()
             self.io.run(self.phys, self.forces, ii, outertime)
             self.phys.updateCOM_Momenta()
             #self.phys.app.energies.clear()
             #self.forces.forcevec.zero()
       else: # Object
	  setPropagator(self, self.phys, self.forces, propFactory.applyModifiers(propFactory.create(1, outerscheme, outertime, outerforcefield, *chain), outerscheme))
          shake = False
          if (params.has_key('shake') and params['shake'] == 'on'):
              shake = True
              shakeMod = ModifierShake.ModifierShake(0.000001, 30)
              #shakeMod = self.myPropagator.createShakeModifier(0.000001, 30)
              self.myPropagator.adoptPostDriftOrNextModifier(shakeMod)
          rattle = False
          if (params.has_key('rattle') and params['rattle'] == 'on'):
              rattle = True
              rattleMod = self.myPropagator.createRattleModifier(0.02, 30)
              self.myPropagator.adoptPostDriftOrNextModifier(rattleMod)
          executePropagator(self, self.phys, self.forces, self.io, steps)  # Runs the propagator for a number of steps
          if (shake):
             self.myPropagator.removeModifier(shakeMod)
          if (rattle):
             self.myPropagator.removeModifier(rattleMod)
       self.phys.updateCOM_Momenta()
           
