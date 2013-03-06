import TopologyUtilities
from ScalarStructure import *

######################################################################
# Class Energies
# Function: 
#####################################################################

class Energies(ScalarStructure):
   """
   Holds the system energies
   These include: potential (bond, angle, ...) and kinetic.

   Inherits from precompiled class ScalarStructure
   This manages energies in a table.  The inherited method getTable()
   will index this table.
   
   This provides several inherited methods which are used:
     - potentialEnergy()
     - intoAdd()
     - intoSubtract()
     - molecularVirial()
     - computeVirial()
   
   """
   def initialize(self, phys):
      self.phys = phys #: Physical object
      
   def computeMolecularVirial(self):
      """
      Tells the energies structure to include the molecular
      virial tensor when calculating terms.

      @rtype: Energies
      @return: The previous state without the virial.
      """

      # Invoke ScalarStructure's molecularVirial()
      return self.phys.app.energies.molecularVirial(1)


   def computeVirial(self):
      """
      Tells the energies structure to include the 
      virial tensor when calculating terms.

      @rtype: Energies
      @return: The previous state without the virial.
      """

      # Invoke ScalarStructure's virial()
      return self.phys.app.energies.virial(1)


   def addBondEnergy(self, r):
      """
      Accumulate into the bond energy
      
      @type r: float
      @param r: Quantity to accumulate.
      
      """
      self.phys.app.energies.setTable(2, self.bondEnergy()+r)

   def addAngleEnergy(self, r):
      """            
      Accumulate into the angle energy

      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(3, self.angleEnergy()+r)
      
   def addDihedralEnergy(self, r):
      """
      Accumulate into the dihedral energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(4, self.dihedralEnergy()+r)

   def addImproperEnergy(self, r):
      """
      Accumulate into the improper energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(5, self.improperEnergy()+r)

   def addShadowEnergy(self, r):
      """
      Accumulate into the shadow energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(34, self.shadowEnergy()+r)

   def addCoulombEnergy(self, r):
      """
      Accumulate into the electrostatic energy
      
      @type r: float
      @param r: Quantity to accumulate.
      """
      self.phys.app.energies.setTable(0, self.coulombEnergy()+r)
      #self.setTable(1, 4)

   def addLJEnergy(self, r):
      """
      Accumulate into the van der Waals energy

      @type r: float
      @param r: Quantity to accumulate.      
      """
      self.phys.app.energies.setTable(1, self.ljEnergy()+r)

            
   def coulombEnergy(self, phys):
      """
      @rtype: float
      @return: Electrostatic energy
      """

      # Coulomb energy is table index 1
      return phys.app.energies.getTable(0)
   
   def ljEnergy(self, phys):
      """
      @rtype: float
      @return: van der Waals energy
      """

      # LJ energy is table index 2
      return phys.app.energies.getTable(1)

   def bondEnergy(self, phys):
      """
      @rtype: float
      @return: Energy due to two-atom bond deviations from equilibrium lengths.
      """

      # Bond energy is table index 3
      return phys.app.energies.getTable(2)

   def angleEnergy(self, phys):
      """
      @rtype: float
      @return: Energy due to three-atom angle deviations from equilibrium angles.
      """

      # Angle energy is table index 4
      return phys.app.energies.getTable(3)

   def dihedralEnergy(self, phys):
      """
      @rtype: float
      @return: Energy due to four-atom dihedral deviations from equilibrium angles.
      """

      # Dihedral energy is table index 5
      return phys.app.energies.getTable(4)
   
   def improperEnergy(self, phys):
      """
      @rtype: float
      @return: Energy due to four-atom improper deviations from equilibrium angles.
      """

      # Improper energy is table index 6
      return phys.app.energies.getTable(5)

   def shadowEnergy(self, phys):
      """
      @rtype: float
      @return: Shadow energy.
      """

      # Shadow energy is table index 4
      return phys.app.energies.getTable(34)

   def kineticEnergy(self, phys):
      """
      @type phys: Physical
      @param phys: Physical system.
      
      @rtype: float
      @return: Kinetic energy, as a sum of 0.5*m*v^2 for each atom.
      """
      return TopologyUtilities.kineticEnergy(phys.myTop, phys.velvec)

   def potentialEnergy(self, phys):
      """
      @type phys: Physical
      @param phys: Physical system.
      
      @rtype: float
      @return: Potential energy
      """
      return phys.app.energies.potentialEnergy()

   def totalEnergy(self, phys):
      """
      @type phys: Physical
      @param phys: Physical system.
      
      @rtype: float
      @return: Total energy, as a sum of potential and kinetic.
      """

      # potentialEnergy() is a member function of ScalarStructure
      return phys.app.energies.potentialEnergy()+TopologyUtilities.kineticEnergy(phys.myTop, phys.velvec)
   
