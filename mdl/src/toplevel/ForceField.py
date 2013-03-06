from ForceGroup import *
import sys
import PySystemForce
class ForceField(ForceGroup):
   """
   A holder for the following data:
     1. Set of forces to evaluate.
     2. Algorithms to use for evaluation.
     3. Switching function(s) to apply, if any.
     4. Extra parameters, depending on the particular forces.
   """
   # DEFAULT DESTRUCTOR WILL DELETE FORCE ARRAY
   # THIS DESTRUCTION IS HANDLED IN THE BACK END
   # SO WE DON'T DO IT AGAIN
   #def __del__(self):
   #   """
   #   Destructor.
   #   """
   #   print "DEF DEL"
   #   pass

   def setDefaults(self):
      """
      Set default values for all parameters, for all force evaluation algorithms.  This includes bonded and nonbonded, pairwise and fast electrostatic.
      """
      self.params = {

       'LennardJones':{'algorithm':'SimpleFull',
                       'switching':'Universal',
                       'cutoff':-1},

       'Coulomb':{'algorithm':'SimpleFull',
                  'switching':'Universal',
                  'cutoff':-1},

       # The 'type' parameter of 'LennardJonesCoulomb' can be set to:
       # 'original', 'EwaldReal', 'EwaldRealTable', 'MagneticDipole',
       # 'MagneticDipoleRealTable'
       'LennardJonesCoulomb':{'algorithm':'SimpleFull',
                              'switching':['Universal','Universal'],
                              'cutoff':-1,
                              'type':'original'},

       'CoulombDiElec':{'algorithm':'SimpleFull',
                        'switching':'Universal',
                        'cutoff':-1},

       # Harmonic Dihedral forces are an exception, in that
       # there can be more than one per force evaluation.
       # i.e., the user may wish to constrain different dihedrals.
       # Thus each parameter is a list of values.
       'HarmonicDihedral':{'kbias':[],          # scaling factor
                           'dihedralnum':[],    # dihedral index
                           'angle':[]}          # constraint angle in radians
       } #: Mapping from parameter names to default values for bonded and nonbonded forces, except fast electrostatics


      # A dirty bit.  If this is set, on a call to propagate()
      # the forces will be rebuilt with the above parameters.
      # Upon each rebuild this bit is cleared, upon each changing
      # of any of the above parameters the bit is set.
      self.dirty = 1  #: Dirty bit, allows for lazy building of force objects upon system propagation
      self.gbsa = False #: Flag telling whether or not to do generalized Born
         
      ##################################################
      # USER-ACCESSIBLE STRUCTURES

      ###################################################################
      # ARRAY OF FORCE TYPES.
      # CAN CONTAIN THE FOLLOWING VALUES:
      # b - BOND, a - ANGLE, d - DIHEDRAL, i - IMPROPER
      # l - LENNARDJONES, c - COULOMB
      # lc - LENNARDJONES/COULOMB TOGETHER PAIRWISE EVALUATION (Default)
      self.forcetypes = [] #: List of force types, can contain 'b' (bond), 'a' (angle), 'd' (dihedral), 'i' (improper), 'h' (harmonic dihedral), 'l' (van der Waals), 'c' (electrostatic), 'e' (implicit solvation dielectric scaling), 'lc' (coupled vdW and electrostatic)
      self.pythonforces = [] #: List of Python-prototyped orce objects
      ################################################################### 
      
      ###################################################################
      # DEFAULT WILL CREATE AN EMPTY FORCE FIELD, BUT IF THE CHARMM
      # FLAG IS SET, WE MUST ADD THE SIX CORE FORCES
      if (self.charmm):
         self.bondedForces("badi")
         self.nonbondedForces("lc")
      ###################################################################   

   def __setattr__(self, s, t):
      if (s == 'params'):   # IF WE CHANGE params or fastelectro, it's dirty
         self.dirty = 1
      self.__dict__[s] = t
      
   # REMOVE ALL BONDED FORCES FROM THE EVALUATION
   def removeBondedForces(self):
      """
      Remove all bonded (bond, angle, dihedral, improper, harmonic dihedral)
      forces from the force field.
      """
      if (self.findForce('b') != -1): self.forcetypes.remove('b')
      if (self.findForce('a') != -1): self.forcetypes.remove('a')
      if (self.findForce('d') != -1): self.forcetypes.remove('d')
      if (self.findForce('i') != -1): self.forcetypes.remove('i')
      if (self.findForce('h') != -1):  # Can be multiple harmonic dihedral forces
         while (self.forcetypes.count('h') != 0): self.forcetypes.remove('h')


   # REMOVE ALL NONBONDED FORCES FROM THE EVALUATION
   def removeNonbondedForces(self):
      """
      Remove all nonbonded forces (van der Waals, electrostatic, magnetic dipole) from the force field.
      """
      if (self.findForce('l') != -1): self.forcetypes.remove('l')
      if (self.findForce('c') != -1): self.forcetypes.remove('c')
      if (self.findForce('lc') != -1): self.forcetypes.remove('lc')
      if (self.findForce('e') != -1): self.forcetypes.remove('e')

      
   # ADD BONDED FORCES ACCORDING TO THE PASSED STRING
   # b=BOND, a=ANGLE, d=DIHEDRAL, i=IMPROPER
   # EXAMPLE INVOCATION: bondedForces("ba")
   def bondedForces(self, inputstring):
      """
      Add bonded forces contained in the input string ('b', 'a', 'd', 'i', or 'h')

      @type inputstring: string
      @param inputstring: Contains the characters representing the bonded forces to instantiate
      """
      self.removeBondedForces()
      if (inputstring.find('b') != -1): self.forcetypes.append('b')
      if (inputstring.find('a') != -1): self.forcetypes.append('a')
      if (inputstring.find('d') != -1): self.forcetypes.append('d')
      if (inputstring.find('i') != -1): self.forcetypes.append('i')
      if (inputstring.count('h') != 0):
          for i in range(0, inputstring.count('h')): self.forcetypes.append('h')
      self.dirty = 1

   # ADD NONBONDED FORCES ACCORDING TO THE PASSED STRING
   # l=LENNARDJONES, c=COULOMB, m=MAGNETIC DIPOLE
   # EXAMPLE INVOCATION: nonbondedForces("lc")
   def nonbondedForces(self, inputstring):
      """
      Add nonbonded forces contained in the input string ('l', 'c', 'm', or 'e').  If 'l' and 'c' are included a coupled vdW-electrostatic is instantiated for maximum performance.

      @type inputstring: string
      @param inputstring: Contains the characters representing the nonbonded forces to instantiate
      """
      self.removeNonbondedForces()
      if (inputstring.find('l') != -1 and
          inputstring.find('c') != -1): self.forcetypes.append('lc')
      else:
          if (inputstring.find('l') != -1): self.forcetypes.append('l')
          if (inputstring.find('c') != -1): self.forcetypes.append('c')
      if (inputstring.find('e') != -1): self.forcetypes.append('e')
      self.dirty = 1

   # FIND A FORCE AND RETURN ITS INDEX
   # IF NOT FOUND, RETURN -1.
   def findForce(self, t):
      """
      Search for a force in the array of force types.

      @type t: char
      @param t: Type of force ('b', 'a', etc. see above).  Force can be bonded or nonbonded.

      @rtype: int
      @return: Index of this type of force in the types array; -1 if not found.
      """
      for ii in range(0, len(self.forcetypes)):
          if (self.forcetypes[ii] == t): return ii
      return -1


   def addPythonForce(self, pyforce):
      """
      Add a Python-prototyped force object for evaluation.

      @type pyforce: PySystemForce
      @param pyforce: An instance of the Python-prototyped force.
      """
      self.pythonforces.append(pyforce)


   # BREAK THE LENNARDJONES/COULOMB EVALUATION OF THE PAIRS
   # TO SEPARATE COMPUTATIONS.  IF THEY ARE ALREADY SEPARATE
   # DO NOTHING
   def breakLennardJonesCoulombForce(self):
      """
      Breaks a coupled van der Waals - electrostatic force into individual calculations.  This is invoked if the user specifies a different algorithm for van der Waals and electrostatic - for example if one is direct and the other uses a cutoff; they will not use the same set of atom pairs.
      """
      pos = self.findForce('lc')
      if (pos != -1):
         self.forcetypes.remove('lc')
         self.forcetypes.append('l')
         self.forcetypes.append('c')

   # THESE LAST TWO FUNCTIONS
   # CAN EVALUATE FORCES DIRECTLY THROUGH THE FORCEFIELD OBJECT
   # ONLY NEED TO USE THEM WHEN CONSTRUCTING A PROPAGATOR FUNCTION
   # (i.e. NOT A CLASS WHICH INHERITS FROM STS OR MTS).  OTHERWISE,
   # DO NOT EVEN USE THEM.
   def calculateForces(self, phys, forces):
      """
      Calculate all forces (except molly) in the force field, and populate the passed MDL forces object.  This should be called from a Python-prototyped propagation function.

      @type phys: Physical
      @param phys: The Physical system.

      @type forces: Forces
      @param forces: MDL Forces object.
      """
      if (forces.forcevec.size() != phys.numAtoms()):
         forces.forcevec.resize(phys.numAtoms())
      forces.forcevec.zero()
      #print len(phys.positions)
      #print phys.app.positions.size()
      phys.app.energies.clear()
      self.phys = forces.phys = phys
      #phys.posvec.setC(phys.positions)
      self.evaluateSystemForces(phys.app, forces.forcevec)
      #sys.exit(1)
      self.evaluateExtendedForces(phys.app, forces.forcevec)

   def build(self):
    """
    Using an MDL force factory, instantiate all force objects stored in the forcetypes data member of each force field.  These will be SWIG-wrapped objects and are appended to the forcearray data member of each force field.
    """
    # TMC 1-13-08: I think this can be improved.  Future:
    # 1. Have a build() member of the ForceField class.
    # 2. Make the ForceFactory a singleton.
    # 3. Give the ForceFactory a method which maps force characters to creation functions.  In this way there are multiple mapping levels.
 
    self.forceFactory.hd = 0

    if (self.params['LennardJones'] !=
        self.params['Coulomb']):
            self.breakLennardJonesCoulombForce()
    self.dirty = 0

    self.forcearray = []

    bornflag=-1
    for forcetype in self.forcetypes:
          if (forcetype == 'b'):
              self.forcearray.append(self.forceFactory.createBondForce(self.bc))
          elif (forcetype == 'a'):
              self.forcearray.append(self.forceFactory.createAngleForce(self.bc))
          elif (forcetype == 'd'):
              self.forcearray.append(self.forceFactory.createDihedralForce(self.bc))
          elif (forcetype == 'i'):
              self.forcearray.append(self.forceFactory.createImproperForce(self.bc))
          elif (forcetype == 'l'):
              self.forcearray.append(self.forceFactory.createLennardJonesForce(self.bc, self.params['LennardJones']))
          elif (forcetype == 'c'):
              if (self.params['Coulomb']['algorithm'] == 'SCPISM'):
                 self.phys.myTop.doSCPISM = 1
                 self.phys.build()
                 #if (not self.params['Coulomb'].has_key('NoBorn')):
                 #   print "CREATING BORN FORCE"
                 #   self.forcearray.append(self.forceFactory.createBornForce(self.bc, self.params['Coulomb']))
                 #   bornflag = len(self.forcearray)-1
                 #else:
                 #    print "NOT CREATING BORN FORCE"
              if (self.params['Coulomb']['algorithm'] == 'GB' or
	          self.params['Coulomb']['algorithm'] == 'GBACE'):
		  self.forcearray.append(self.forceFactory.createBornBurial())
                  self.addForce(self.forcearray[len(self.forcearray)-1])
	          self.forcearray.append(self.forceFactory.createBornRadii())
                  self.addForce(self.forcearray[len(self.forcearray)-1])
              if (not self.params['Coulomb'].has_key('OnlyBorn')):      
                 self.forcearray.append(self.forceFactory.createCoulombForce(self.bc, self.params['Coulomb']))
          elif (forcetype == 'e'):
              self.forcearray.append(self.forceFactory.createCoulombDiElecForce(self.bc, self.params['CoulombDiElec']))
          elif (forcetype == 'lc'):
             self.forcearray.append(self.forceFactory.createLennardJonesCoulombForce(self.bc, self.params['LennardJonesCoulomb']))
          elif (forcetype == 'h'):
              self.forcearray.append(self.forceFactory.createHarmDihedralForce(self.bc, self.params['HarmonicDihedral']))
          self.addForce(self.forcearray[len(self.forcearray)-1])
	  import ForceGroup
	  ForceGroup._swig_setattr_nondynamic(self.forcearray[len(self.forcearray)-1], ForceGroup.Force, "thisown", 0)

    if (bornflag != -1): self.forcetypes.insert(bornflag, 'c')
    for pyforce in self.pythonforces:
         self.forcearray.append(PySystemForce.PySystemForce(pyforce))
	 PySystemForce._swig_setattr_nondynamic(self.forcearray[len(self.forcearray)-1], PySystemForce.PySystemForce, "thisown", 0)
         self.addSystemForce(self.forcearray[len(self.forcearray)-1])
