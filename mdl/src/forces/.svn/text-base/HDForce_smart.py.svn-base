import numpy
import umath
from Potential import *
import sys

def norm(v3d):
    return numpy.sqrt(v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2])

def norm2(v3d):
    return v3d[0]*v3d[0]+v3d[1]*v3d[1]+v3d[2]*v3d[2]
  
class HDForce:
  def __init__(self, phys, forces, phi, num, k):
      self.phys = phys
      self.forces = forces
      self.phi = phi
      self.dihedral = num
      self.k = k

  def dihed(self, x):
      atomI = self.phys.dihedral(self.dihedral-1).atom1 - 1
      atomJ = self.phys.dihedral(self.dihedral-1).atom2 - 1
      atomK = self.phys.dihedral(self.dihedral-1).atom3 - 1
      atomL = self.phys.dihedral(self.dihedral-1).atom4 - 1
      r12 = self.phys.positions[atomJ*3:atomJ*3+3] - self.phys.positions[atomI*3:atomI*3+3]
      r23 = self.phys.positions[atomK*3:atomK*3+3] - self.phys.positions[atomJ*3:atomJ*3+3]
      r34 = self.phys.positions[atomL*3:atomL*3+3] - self.phys.positions[atomK*3:atomK*3+3]
      a = numpy.cross(r12, r23)
      b = numpy.cross(r23, r34)
      c = numpy.cross(r23, a)
      cosPhi = numpy.dot(a,b) / (norm(a)*norm(b))
      sinPhi = numpy.dot(c,b) / (norm(c)*norm(b))
      return -math.atan(sinPhi/cosPhi)
      
  def V(self, x):
      dd = self.dihed(x) - self.phi
      if (dd < -numpy.pi):
        dd += 2 * numpy.pi
      elif (dd > numpy.pi):
        dd -= 2 * numpy.pi
      return self.k*dd**2

  def evaluate(self):
      atomI = self.phys.dihedral(self.dihedral-1).atom1 - 1
      atomJ = self.phys.dihedral(self.dihedral-1).atom2 - 1
      atomK = self.phys.dihedral(self.dihedral-1).atom3 - 1
      atomL = self.phys.dihedral(self.dihedral-1).atom4 - 1
      x1, x2, x3, x4 = DerivVectors(Vector(*tuple(self.phys.positions[atomI*3:atomI*3+3])),
                                    Vector(*tuple(self.phys.positions[atomJ*3:atomJ*3+3])),
                                    Vector(*tuple(self.phys.positions[atomK*3:atomK*3+3])),
                                    Vector(*tuple(self.phys.positions[atomL*3:atomL*3+3])))

      a = Vector.cross(x1-x2, x2-x3)
      b = Vector.cross(x2-x3, x3-x4)
      c = Vector.cross(x2-x3, a)
      cosPhi = a*b / (a.length()*b.length())
      sinPhi = c*b / (c.length()*b.length())
      dihed = -(sinPhi.arctan2(cosPhi))
      dd = dihed - self.phi
      if (dd < -numpy.pi):
        dd += 2 * numpy.pi
      elif (dd > numpy.pi):
        dd -= 2 * numpy.pi
      e = self.k*(dd**2)

      energy, forces = EnergyGradients(e, 4)
      self.forces.energies.addDihedralEnergy(energy)
      self.forces.force[atomI*3:atomI*3+3] -= forces[0].data
      self.forces.force[atomJ*3:atomJ*3+3] -= forces[1].data
      self.forces.force[atomK*3:atomK*3+3] -= forces[2].data
      self.forces.force[atomL*3:atomL*3+3] -= forces[3].data
      sys.stdout.flush()

