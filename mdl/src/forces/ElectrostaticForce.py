import numpy
import math
#from Potential import *
import sys
import Constants
def norm(v3d):
    return numpy.sqrt(v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2])

def norm2(v3d):
    return v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2]
  
class ElectrostaticForce:
  """
  Implement a harmonic dihedral constraining potential:
  U(x) = k*(phi - phi0)^2
  """
  
  def __init__(self, phys, forces):
    """
    Initialize an object of type HDForce
    
    @type phys: Physical
    @param phys: The physical system.

    @type forces: Forces
    @param forces: MDL Forces object

    """
    self.phys = phys
    self.forces = forces

    n = self.phys.numAtoms()
    self.q = []
    
    for i in range (0, n):
        self.q.append(self.phys.charge(i+1))
      
  def eval(self):
    """
    Modify energy and force vector to include this force term.
    """

    n = self.phys.numAtoms()
    for i in range (0, n):
        for j in range (i+1, n):
            rij = self.phys.positions[j*3:j*3+3] - self.phys.positions[i*3:i*3+3]
            rdistsquared = 1.0 / norm2(rij)
            en = self.q[i]*self.q[j]*numpy.sqrt(rdistsquared)
            self.forces.energies.addCoulombEnergy(en)
            fo = en*rdistsquared
            self.forces.force[i:i+3] -= fo
            self.forces.force[j:j+3] += fo

