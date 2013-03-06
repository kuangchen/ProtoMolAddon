from Energies import *
import Vector3DBlock
from ForceField import *
from ForceFactory import *
import MathUtilities
import numpy
from PySystemForce import *

class Forces:
   """
   Contains the atomic force vector, and a structure
   to hold system energies.
   """
   def __init__(self):
      #####################################################################
      # USER-ACCESSIBLE STRUCTURES
      self.myForceFields = list()   #: Array of MDL force fields
      self.energies = Energies()    #: Holds system energies
      # NOTE: self.forces is a numpy array which holds system forces.
      #####################################################################

      #########################################################
      # NOTE: USER SHOULD NOT TOUCH THIS!
      # THIS IS A WRAPPER FOR THE FORCE ARRAY
      # WHICH IS ACCESSIBLE THROUGH self.force
      self.__dict__['forcevec'] = Vector3DBlock.Vector3DBlock()
      #self.force = 0  #: Atomic force vector.
      #########################################################
 
      self.accept = 0
      self.reject = 0
         
   # SPECIAL ACCESSOR FOR self.forces
   # TO GET DATA FROM WRAPPERS   
   def __getattr__(self, name):
      if (name == 'force'):
         return self.__dict__['forcevec'].getC()
      else:
         return self.__dict__[name]

   # SPECIAL ASSIGNMENT FOR self.forces
   # TO SET DATA IN WRAPPERS 
   def __setattr__(self, name, val):
      if (name == 'force' and self.__dict__.has_key('force')):
         self.__dict__['forcevec'].setC(val)
      else:
         self.__dict__[name] = val

   def __str__(self):
      return "Forces.Forces"
   
   def __repr__(self):
      return "Forces.Forces"


   def initializeEnergies(self, app):
      """
      Initialize the energies structure.
      """
      self.energies = app.energies
      
   # RESET SYSTEM STATE TO DEFAULTS
   def reset(self):
      """
      Reset data members to default values.
      """
      self.myForceFields = list()
      self.energies = Energies()
      self.forces.fill(0)

   # CREATE A FORCE FIELD AND APPEND IT TO THE END
   # OF THE ARRAY.
   # ALSO RETURN THE NEW FORCE FIELD POINTER.
   def makeForceField(self, phys, *args):
      """
      Create a new MDL force field.

      @type phys: Physical
      @param phys: The physical system.

      @type args: tuple
      @param args: Python tuple, may have no values or a string \"charmm\"

      @rtype: ForceField
      @return: Newly instantiated MDL force field.
      """
      ff = ForceField()
      ff.charmm = 0
      ff.forceFactory = ForceFactory()
      if (args.__len__() > 0):
         ff.charmm = 1
      ff.theforces = self
      self.phys = ff.phys = phys
      ff.bc = phys.bc
      ff.setDefaults()
      self.energies.initialize(phys)
      ff.thisown = 0
      #self.myForceFields.append(ff)
      return ff
      #return self.myForceFields[self.myForceFields.__len__()-1]


   # NOTE: THIS WILL PERMANENTLY DELETE THE PASSED FORCE FIELD OBJECT
   # FROM THE ARRAY!  USE WITH CAUTION.
   def removeForceField(self, ff):
      """
      Remove the passed force field from the member list.

      @type ff: ForceField
      @param ff: MDL force field.      
      """
      for i in range(0, self.myForceFields.__len__()):
         if (self.myForceFields[i] == ff):
	    self.myForceFields.remove(ff)
            break

   def randomForce(self, phys, seed):
        """
        Compute a random (x, y, z) force.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type seed: integer
        @param seed: Random number generator seed.

        @rtype: numpy.ndarray
        @return: The random force as a three-element array (x,y,z)
        """
        retval = numpy.ndarray(phys.numAtoms()*3)
        #tmp = numpy.ndarray(3)

        ff = 0        
        while (ff < phys.numAtoms()*3):
           retval[ff] = MathUtilities.randomGaussianNumber(seed)
           retval[ff+1] = MathUtilities.randomGaussianNumber(seed)
           retval[ff+2] = MathUtilities.randomGaussianNumber(seed)
           ff += 3
        return retval

