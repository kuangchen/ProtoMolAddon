import numpy
import math
import sys
from PySystemForce import *

def norm(v3d):
    return numpy.sqrt(v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2])

def norm2(v3d):
    return v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2]
  
class HDForce:
  """
  Implement a harmonic dihedral constraining potential:
  U(x) = k*(phi - phi0)^2
  """
  
  def __init__(self, phys, forces, phi, num, k):
    """
    Initialize an object of type HDForce
    
    @type phys: Physical
    @param phys: The physical system.

    @type forces: Forces
    @param forces: MDL Forces object
    
    @type phi: float
    @param phi: Target dihedral value in radians.

    @type num: integer
    @param num: Dihedral number to constrain.

    @type k: float
    @param k: Scaling factor for constrain.
    """
    self.phys = phys
    self.forces = forces
    self.phi = phi
    self.dihedral = num
    self.k = k
      
  def eval(self):
    """
    Modify energy and force vector to include this force term.
    """
    # COMPUTE DIHEDRAL DIFFERENCE
    diff = self.phys.dihedral(self.dihedral) - self.phi;
    if (diff < -numpy.pi):
      diff += 2 * numpy.pi;
    elif (diff > numpy.pi):
      diff -= 2 * numpy.pi;

    # ADD ENERGY
    self.forces.energies.addDihedralEnergy(self.k * diff * diff)
    # COMPUTE FORCES
    atomI = self.phys.dihedral(self.dihedral).atom1 - 1
    atomJ = self.phys.dihedral(self.dihedral).atom2 - 1
    atomK = self.phys.dihedral(self.dihedral).atom3 - 1
    atomL = self.phys.dihedral(self.dihedral).atom4 - 1
    rij = self.phys.positions[atomJ*3:atomJ*3+3] - self.phys.positions[atomI*3:atomI*3+3]
    rkj = self.phys.positions[atomJ*3:atomJ*3+3] - self.phys.positions[atomK*3:atomK*3+3]
    rkl = self.phys.positions[atomL*3:atomL*3+3] - self.phys.positions[atomK*3:atomK*3+3]
    m = numpy.cross(rij, rkj)
    n = numpy.cross(rkj, rkl)
    dVdPhi = 2 * self.k * diff

    fi = m * (-dVdPhi * norm(rkj) / norm2(m));
    fl = n * (dVdPhi * norm(rkj) / norm2(n));
    fj = fi * (-1 + numpy.dot(rij, rkj)/norm2(rkj)) - fl * (numpy.dot(rkl, rkj)/norm2(rkj));
    fk = - (fi + fj + fl);

    self.forces.force[atomI*3:atomI*3+3] += fi
    self.forces.force[atomJ*3:atomJ*3+3] += fj
    self.forces.force[atomK*3:atomK*3+3] += fk
    self.forces.force[atomL*3:atomL*3+3] += fl

