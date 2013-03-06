import warnings
warnings.filterwarnings(action='ignore',
                        message='.*has API version.*',
                        category=RuntimeWarning)

import Forces
import Constants
import Vector3DBlock
import PARReader
import PSFReader
import PDBReader
import XYZReader
import EigenvectorReader
import MathUtilities
import TopologyUtilities
import GenericTopology
#import _TopologyUtilities
import _GenericTopology
import numpy
import sys
import ProtoMolApp
def deepcopy(a, b):
   if (len(b) != len(a)):
      b.resize(len(a))
   for i in range (0, len(a)):
      b[i] = a[i]


class MDVec(numpy.ndarray):
   def __eq__(self, b):
    if (len(self) != len(b)):
      self.resize(len(b))
    for i in range (0, len(b)):
      self[i] = b[i]
      

class Atom:
   """
   An atom in the system.
   """
   def __init__(self, n, s, rs, rn, an, at, c, m):
      self.number = n   #: Atom number
      self.seg_id = s   #: Segment identifier
      self.residue_sequence = rs  #: Residue sequence
      self.residue_name = rn     #: Residue name
      self.atom_name = an        #: Atom name
      self.atom_type = at        #: Atom type
      self.charge = c            #: Charge [e]
      self.mass = m              #: Mass [amu]

class Bond:
   """
   A two-atom bond.
   """
   def __init__(self, n, a1, a2):
      self.number = n    #: Bond number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2

class Angle:
   """
   A three-atom angle.
   """
   def __init__(self, n, a1, a2, a3):
      self.number = n    #: Bond number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2
      self.atom3 = a3    #: Atom index 3

class Dihedral:
   """
   A four-atom dihedral.
   """
   def __init__(self, n, a1, a2, a3, a4):
      self.number = n    #: Dihedral number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2
      self.atom3 = a3    #: Atom index 3
      self.atom4 = a4    #: Atom index 4

class Improper:
   """
   A four-atom improper.
   """
   def __init__(self, n, a1, a2, a3, a4):
      self.number = n    #: Bond number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2
      self.atom3 = a3    #: Atom index 3
      self.atom4 = a4    #: Atom index 4

class HDonor:
   """
   An H+ donor
   """
   def __init__(self, n, a1, a2):
      self.number = n    #: Donor number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2

class HAcceptor:
   """
   An H+ acceptor
   """
   def __init__(self, n, a1, a2):
      self.number = n    #: Acceptor number
      self.atom1 = a1    #: Atom index 1
      self.atom2 = a2    #: Atom index 2


global apps
apps = []     

class Physical:
   """
   Defines a physical system: positions, velocities, temperature,
   boundary conditions, etc.
   """


   def __init__(self):
      #####################################################################
      # USER-ACCESSIBLE STRUCTURES
      self.seed = 1234               #: Random number seed
      self.exclude = "1-4"           #: Exclusion Pairs
      self.cellsize = 6.5            #: Cell size
      self.bc = "Periodic"           #: Boundary conditions
      self.remcom = 0                #: Remove COM motion?
      self.remang = 0                #: Remove angular momentum?
      self.defaultCBV = True         #: Use default cell basis vectors (PBC)
      self.time = 0                  #: Current time
      self.cB1 = numpy.ndarray(3)    #: Cell basis vector 1
      self.cB1.fill(0)
      self.cB2 = numpy.ndarray(3)    #: Cell basis vector 2
      self.cB2.fill(0)
      self.cB3 = numpy.ndarray(3)    #: Cell basis vector 3
      self.cB3.fill(0)
      self.cO = numpy.ndarray(3)     #: Cell origin
      self.cO.fill(0)
      self.temperature = 300         #: Kelvin temperature
      self.masses = numpy.ndarray(0)      #: Diagonal mass matrix
      self.invmasses = numpy.ndarray(0)   #: Diagonal inverse mass matrix
      self.masssum = 0                    #: Sum over all atomic masses
      self.accept = 0
      self.reject = 0
      ProtoMolApp.ProtoMolApp.turnOffHints()
      # NOTE: self.positions and self.velocities are also
      # available numpy arrays.
      #####################################################################

      ############################################
      # NOTE: USER SHOULD NOT TOUCH THESE!
      # THESE ARE NUMPY ARRAY WRAPPERS
      # USER SHOULD ACCESS THROUGH
      # self.positions and self.velocities
      self.__dict__['myTop'] = GenericTopology.T_Periodic()
      self.__dict__['posvec'] = Vector3DBlock.Vector3DBlock()
      self.__dict__['velvec'] = Vector3DBlock.Vector3DBlock()
      self.__dict__['myPAR'] = PARReader.PAR()
      self.__dict__['myPSF'] = PSFReader.PSF()
      self.__dict__['myPDB'] = PDBReader.PDB()
      self.__dict__['myXYZ'] = XYZReader.XYZ()
      self.__dict__['myEig'] = EigenvectorReader.EigenvectorInfo()
      ############################################

      #self.positions = numpy.ndarray(0)   #: Atomic position vector
      #self.velocities = numpy.ndarray(0)  #: Atomic velocity vector

      self.dirty = 1   #: Dirty bit
   
   def __del__(self):
      pass

   def __str__(self):
      return "Physical.Physical"

   def __repr__(self):
      return "Physical.Physical"

   # Copy which avoids object assignment
   def copy(self, retval, forces="", dt=-1):
      """
      Perform a deep copy, avoid reference assignment
      """
      #retval = Physical()
      retval.seed = self.seed
      retval.exclude = self.exclude
      retval.cellsize = self.cellsize
      retval.bc = self.bc
      retval.defaultCBV = self.defaultCBV
      retval.remcom = -1
      retval.remang = -1
      retval.defaultCBV = self.defaultCBV
      retval.myTop = self.myTop
      #if (retval.bc == "Periodic"):
      #   retval.myTop = GenericTopology.T_Periodic()
      #else:
      #   retval.myTop = GenericTopology.T_Vacuum()
      retval.cB1 = self.cB1.copy()
      retval.cB2 = self.cB2.copy()
      retval.cB3 = self.cB3.copy()
      retval.cO = self.cO.copy()
      retval.temperature = self.temperature
      retval.masses = self.masses.copy()
      retval.invmasses = self.invmasses.copy()
      retval.masssum = self.masssum
      retval.posvec.resize(self.posvec.size())
      for i in range(len(self.positions)):
         retval.positions[i] = self.positions[i]
      retval.velvec.resize(self.velvec.size())
      for i in range(len(self.velocities)):
         retval.velocities[i] = self.velocities[i]
      retval.myPAR = self.myPAR
      retval.myPSF = self.myPSF
      retval.myPDB = self.myPDB
      retval.myEig = self.myEig
      retval.accept = self.accept
      retval.reject = self.reject

      if (dt != -1):
         retval.app = ProtoMolApp.ProtoMolApp()
         f = Forces.Forces()
      	 retval.app.makeApp(retval.myTop, retval.posvec, retval.velvec, f.energies, dt)
         apps.append(retval.app)
      #retval.build()
      #return retval
      

   def deepcopy(self, forces, dt):
         retval = self.copy()
         retval.app = ProtoMolApp.ProtoMolApp()
         retval.app.makeApp(retval.myTop, retval.posvec, retval.velvec, forces.energies, dt)

      
   # SPECIAL ACCESSOR FOR self.positions or self.velocities
   # TO GET DATA FROM WRAPPERS
   def __getattr__(self, name):
      if (name == 'positions'):
	 return (self.__dict__['posvec'].getC()).view(MDVec)
      elif (name == 'velocities'):
         return (self.__dict__['velvec'].getC()).view(MDVec)
      elif (name == 'time'):
         return self.myTop.time
      else:
         print name
         return self.__dict__[name]

   # SPECIAL ASSIGNMENT FOR self.positions or self.velocities
   # TO SET DATA IN WRAPPERS
   def __setattr__(self, name, val):
      firsttime = False
      if (not self.__dict__.has_key(name)):
         firsttime = True
      if (name == 'positions'):
	 for i in range (0, len(self.positions)):
	    self.positions[i] = val[i]
	 #self.__dict__['posvec'].setC(val)
      elif (name == 'velocities'):
	 for i in range (0, len(self.velocities)):
	    self.velocities[i] = val[i]
	 #self.__dict__['velvec'].setC(val)
      elif (name == 'bc'):
         self.__dict__['bc'] = val         
         if (val == "Vacuum"):
           self.myTop = GenericTopology.T_Vacuum()
         else:
           self.myTop = GenericTopology.T_Periodic()
         if (not firsttime):
            self.build()
      elif (name == 'time'):
         val /= Constants.invTimeFactor()
         self.myTop.time = val
      elif (name == 'seed' or
            name == 'exclude' or
            name == 'cellsize' or
            name == 'remcom' or
            name == 'remang'):
         self.__dict__[name] = val
         if (not firsttime):
            self.build()
      else:
         self.__dict__[name] = val

   # RESETS SIMULATION STATE TO DEFAULTS
   # ALSO RESETS TIME TO ZERO
   def reset(self):
      """
      Reset all member variables to default values
      """
      self.__dict__['seed'] = 1234               # Random number seed
      self.__dict__['exclude'] = "1-4"           # Exclusion Pairs
      self.__dict__['cellsize'] = 6.5            # Cell size
      self.__dict__['bc'] = "Periodic"           # Boundary conditions
      self.__dict__['remcom'] = 0                # Remove COM motion?
      self.__dict__['remang'] = 0                # Remove angular momentum?
      self.__dict__['defaultCBV'] = True         # Use default cell basis vectors (PBC)
      self.__dict__['myTop'] = Topology.T_Periodic()
      self.time = 0
      self.cB1.fill(0)
      self.cB2.fill(0)
      self.cB3.fill(0)
      self.cO.fill(0)
      self.__dict__['temperature'] = 300         # Kelvin temperature

   # SYSTEM PRESSURE.
   def pressure(self, forces):
      """
      Pressure of the system.

      @type forces: Forces
      @param forces: MDL Forces object

      @rtype: float
      @return: System pressure
      """
      TopologyUtilities.addVelocityVirial(forces.energies, self.myTop, self.velvec)
      return TopologyUtilities.velocityVirial(self.myTop, self.velvec).pressure(self.myTop.getVolume(self.posvec))

   # SYSTEM VOLUME (AA^3)
   def volume(self):
      """
      Volume of the system.

      @rtype: float
      @return: System volume
      """
      return self.myTop.getVolume(self.posvec)


   # SIZE OF THE SYSTEM
   def N(self):
      """
      Number of atoms.

      @rtype: int
      @return: Number of atoms.
      """
      return self.numAtoms()
   
   def numAtoms(self):
      """
      Number of atoms.

      @rtype: int
      @return: Number of atoms.
      """
      return self.myPSF.numAtoms()

   # NUMBER OF TWO-ATOM BONDS.
   def numBonds(self):
      """
      Number of bonds.

      @rtype: int
      @return: Number of bonds
      """
      return self.myPSF.numBonds()

   # NUMBER OF THREE-ATOM ANGLES.
   def numAngles(self):
      """
      Number of angles

      @rtype: int
      @return: Number of angles
      """
      return self.myPSF.numAngles()

   # NUMBER OF FOUR-ATOM DIHEDRALS.
   def numDihedrals(self):
      """
      Number of dihedrals

      @rtype: int
      @return: Number of dihedrals
      """      
      return self.myPSF.numDihedrals()

   # NUMBER OF FOUR-ATOM IMPROPERS.
   def numImpropers(self):
      """
      Number of impropers

      @rtype: int
      @return: Number of impropers
      """
      return self.myPSF.numImpropers()

   # NUMBER OF H-BOND DONORS.
   def numDonors(self):
      """
      Number of hydrogen donors (for H+ bonding)

      @rtype: int
      @return: Number of hydrogen donors
      """
      return self.myPSF.numDonors()

   # NUMBER OF H-BOND ACCEPTORS.
   def numAcceptors(self):
      """
      Number of hydrogen acceptors (for H+ bonding)

      @rtype: int
      @return: Number of hydrogen acceptors
      """
      return self.myPSF.numAcceptors()

   # GET ATOM #(index)
   def atom(self, index):
      """
      Get an atom at the passed index.

      @type index: int
      @param index: Atom index (1 to N)

      @rtype: Atom
      @return: The atom at the passed index.
      """
      return Atom(self.myPSF.getAttributeInt("atom", index-1, "number"),
                  self.myPSF.getAttributeString("atom", index-1, "seg_id"),
                  self.myPSF.getAttributeInt("atom", index-1, "residue_sequence"),
                  self.myPSF.getAttributeString("atom", index-1, "residue_name"),
                  self.myPSF.getAttributeString("atom", index-1, "atom_name"),
                  self.myPSF.getAttributeString("atom", index-1, "atom_type"),
                  self.myPSF.getAttributeReal("atom", index-1, "charge"),
                  self.myPSF.getAttributeReal("atom", index-1, "mass"))


   # GET BOND #(index)
   def getBond(self, index):
      """
      Get a bond at the passed index.

      @type index: int
      @param index: Bond index

      @rtype: Bond
      @return: The bond at the passed index.
      """
      return Bond(self.myPSF.getAttributeInt("bond", index-1, "number"),
                  self.myPSF.getAttributeInt("bond", index-1, "atom1"),
                  self.myPSF.getAttributeInt("bond", index-1, "atom2"))

   # GET ANGLE #(index)
   def getAngle(self, index):
      """
      Get an angle at the passed index.

      @type index: int
      @param index: Angle index

      @rtype: Angle
      @return: The angle at the passed index.
      """
      return Angle(self.myPSF.getAttributeInt("angle", index-1, "number"),
                   self.myPSF.getAttributeInt("angle", index-1, "atom1"),
                   self.myPSF.getAttributeInt("angle", index-1, "atom2"),
                   self.myPSF.getAttributeInt("angle", index-1, "atom3"))

   # GET DIHEDRAL #(index)
   def getDihedral(self, index):
      """
      Get a dihedral at the passed index.

      @type index: int
      @param index: Dihedral index

      @rtype: Dihedral
      @return: The dihedral at the passed index.
      """
      return Dihedral(self.myPSF.getAttributeInt("dihedral", index-1, "number"),
                      self.myPSF.getAttributeInt("dihedral", index-1, "atom1"),
                      self.myPSF.getAttributeInt("dihedral", index-1, "atom2"),
                      self.myPSF.getAttributeInt("dihedral", index-1, "atom3"),
                      self.myPSF.getAttributeInt("dihedral", index-1, "atom4"))

   # GET IMPROPER #(index)
   def getImproper(self, index):
      """
      Get an improper at the passed index.

      @type index: int
      @param index: Improper index

      @rtype: Improper
      @return: The improper at the passed index.
      """
      return Improper(self.myPSF.getAttributeInt("improper", index-1, "number"),
                      self.myPSF.getAttributeInt("improper", index-1, "atom1"),
                      self.myPSF.getAttributeInt("improper", index-1, "atom2"),
                      self.myPSF.getAttributeInt("improper", index-1, "atom3"),
                      self.myPSF.getAttributeInt("improper", index-1, "atom4"))


   # GET DONOR #(index)
   def getDonor(self, index):
      """
      Get an H+ donor at the passed index.

      @type index: int
      @param index: H+ donor index

      @rtype: HDonor
      @return: The H+ donor at the passed index.
      """
      return HDonor(self.myPSF.getAttributeInt("donor", index-1, "number"),
                    self.myPSF.getAttributeInt("donor", index-1, "atom1"),
                    self.myPSF.getAttributeInt("donor", index-1, "atom2"))

   # GET ACCEPTOR #(index)
   def getAcceptor(self, index):
      """
      Get an H+ acceptor at the passed index.

      @type index: int
      @param index: H+ acceptor index

      @rtype: HAcceptor
      @return: The H+ acceptor at the passed index.
      """
      return HAcceptor(self.myPSF.getAttributeInt("acceptor", index-1, "number"),
                       self.myPSF.getAttributeInt("acceptor", index-1, "atom1"),
                       self.myPSF.getAttributeInt("acceptor", index-1, "atom2"))



   # GET THE MASS (AMU) OF A SPECIFIC ATOM
   def mass(self, atom):
      """
      Mass of an atom.

      @type atom: int
      @param atom: Atom index (1 to N)

      @rtype: float
      @return: Atom mass [amu]
      """
      return self.myPSF.getAttributeReal("atom", atom-1, "mass")

   def charge(self, atom):
      """
      Charge of an atom.

      @type atom: int
      @param atom: Atom index (1 to N)

      @rtype: float
      @return: Atomic charge [e]
      """
      return self.myPSF.getAttributeReal("atom", atom-1, "charge")

   # SYSTEM TEMPERATURE.
   def temperature(self):
      """
      System temperature (K)

      @rtype: float
      @return: Kelvin temperature
      """
      return TopologyUtilities.temperature(self.myTop, self.velvec)

   def dihedral(self, index):
      """
      Dihedral angle (rad, -PI to PI) at passed index

      @type index: int
      @param index: Dihedral index
      
      @rtype: float
      @return: Dihedral angle in radians
      """
      myPhi = TopologyUtilities.computePhiDihedral(self.myTop, self.posvec, index-1)
      if (myPhi > numpy.pi):
         myPhi -= 2*numpy.pi
      elif (myPhi < -numpy.pi):
         myPhi += 2*numpy.pi
      return myPhi


   def randomVelocity(self, T):
      """
      Assign random velocities.

      @type T: float
      @param T: Kelvin temperature.
      """
      TopologyUtilities.randomVelocity(T, self.myTop, self.velvec, self.seed)

      
   def updateCOM_Momenta(self):
      """
      Update center of mass and angular momentum
      """
      TopologyUtilities.buildMolecularCenterOfMass(self.posvec,self.myTop)
      TopologyUtilities.buildMolecularMomentum(self.velvec,self.myTop)

   def build(self):
      """
      Build the physical data.
      """
      tm = -1
      if (hasattr(self, "myTop")):
          tm = self.myTop.time
      if (self.bc == "Periodic"):
          if (self.defaultCBV):
             # TEMPORARY STRUCTURES USED TO COMPUTE
             # BOUNDING BOX
             v1 = numpy.ndarray(3)
             v2 = numpy.ndarray(3)
             v3 = numpy.ndarray(3)
             v4 = numpy.ndarray(3)
             v1.fill(sys.maxint)
             v2.fill(-sys.maxint)
             # BOUNDING BOX
             i = 0
             while (i < numpy.size(self.positions)):
                if (self.positions[i] < v1[0]):
                   v1[0] = float(self.positions[i])
                if (self.positions[i] > v2[0]):
                   v2[0] = float(self.positions[i])
                if (self.positions[i+1] < v1[1]):
                   v1[1] = float(self.positions[i+1])
                if (self.positions[i+1] > v2[1]):
                   v2[1] = float(self.positions[i+1])
                if (self.positions[i+2] < v1[2]):
                   v1[2] = float(self.positions[i+2])
                if (self.positions[i+2] > v2[2]):
                   v2[2] = float(self.positions[i+2])
                i += 3
             v4.fill(Constants.periodicBoundaryTolerance()/2.0)
             v1 = v1 - v4
             v2 = v2 + v4
             v3 = v2 - v1
             self.cB1[0] = v3[0]
             self.cB2[1] = v3[1]
             self.cB3[2] = v3[2]
             self.cO = v1 + v3 * 0.5
             self.myTop.setBC(self.cB1[0],self.cB1[1],self.cB1[2],self.cB2[0],self.cB2[1],self.cB2[2],self.cB3[0],self.cB3[1],self.cB3[2],self.cO[0],self.cO[1],self.cO[2])
      self.myTop.setExclusion(self.exclude)
      if (self.myPSF.numAtoms() > 0 and hasattr(self.myPAR, 'readFlag')):
	 GenericTopology.buildTopology(self.myTop, self.myPSF, self.myPAR, 0, self.myTop.makeSCPISM())
	 if (numpy.size(self.velocities) == 0):
            # This function actually returns a value, since it is
            # unused the runtime environment detects a memory leak.
            # disown() removes this concern.
            MathUtilities.randomNumber(self.seed)
            # NOTE TO SELF: THIS WAS COMMENTED OUT FOR PM 3
            # IT MAY EFFECT RANDOM NUMBER CONSISTENCY
            #aaa = MathUtilities.randomNumberFirst(self.seed, 1)
	    TopologyUtilities.randomVelocity(self.temperature, self.myTop, self.velvec, self.seed)
	 if (self.remcom >= 0):
	   TopologyUtilities.removeLinearMomentum(self.velvec, self.myTop).disown()
         if (self.remang >= 0):
	   TopologyUtilities.removeAngularMomentum(self.posvec, self.velvec, self.myTop).disown()
      if (self.bc == "Periodic"):
           self.myTop.setCellSize(self.cellsize)
      # COMPUTE INV MASS MATRIX
      temp = list()
      ii = 0
      while ii < self.numAtoms()*3:
          temp.append(0.0)
          ii += 1
      ii = 0
      self.invmasses.resize(self.numAtoms()*3)
      while ii < self.numAtoms()*3:
          im = 1.0 / self.myPSF.getMass(ii/3)
          self.invmasses[ii] = im
          self.invmasses[ii+1] = im
          self.invmasses[ii+2] = im
          ii += 3

      # COMPUTE MASS MATRIX
      temp = list()
      ii = 0
      while ii < self.numAtoms()*3:
          temp.append(0.0)
          ii += 1
      ii = 0
      self.masses.resize(self.numAtoms()*3)
      while ii < self.numAtoms()*3:
          m = self.myPSF.getMass(ii/3)
          temp[ii] = m
          self.masses[ii] = temp[ii]
          temp[ii] = 0.0
          temp[ii+1] = m
          self.masses[ii+1] = temp[ii+1]
          temp[ii+1] = 0.0
          temp[ii+2] = m
          self.masses[ii+2] = temp[ii+2]
          temp[ii+2] = 0.0
          ii += 3


      # COMPUTE MASS SUM
      self.masssum = 0
      ii = 0
      while ii < self.numAtoms():
         self.masssum += self.myPSF.getMass(ii)
         ii += 1

      # SET SELFICAL TIME
      if (tm != -1):
          self.myTop.time = tm

      self.dirty = 0  # clean now!


