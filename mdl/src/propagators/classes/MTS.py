import PyMTSIntegrator
from ForceGroup import *
import ScalarStructure

class MTS(PyMTSIntegrator.PyMTSIntegrator):
   """
   Parent class for all multiple time-stepping propagators.
   """
   def init(self):
      """
      Initialize the propagator.
      """
      return
   def run(self):
      """
      Run the propagator.
      """
      return
   def finish(self):
      """
      Finalize the propagator.
      """
      return


class MOLLY(MTS):
  """
  Parent class for all MOLLY propagators.
  """
  def setMOLLYForceGroups(self):
    """
    Initialize MOLLY force fields and energies.
    """
    self.myMOLLYForcesBonded = ForceGroup()
    self.myMOLLYForcesHBonded = ForceGroup()
    self.myHBondForces = ForceGroup()
    self.mollyEnergy = ScalarStructure.ScalarStructure()
  def calculateBondedForces(self):
    """
    Evaluate mollified bonded forces.
    """
    self.mollyForces.zero()
    self.myMOLLYForcesBonded.evaluateSystemForces(self.getTopo(), self.mollyPos, self.mollyForces, self.mollyEnergy)
