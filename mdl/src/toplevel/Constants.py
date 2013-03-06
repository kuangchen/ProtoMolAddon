######################################################################
# MDL Constants                                                      #
# Many of these convert to SI units for calculations                 #
######################################################################

def boltzmann():
    """
    @rtype: float
    @return: Boltzmann constant, units kcal/mol K^-1
    """
    return 0.001987191

def epsilon():
    """
    @rtype: float
    @return: 10^-14
    """
    return 1*pow(10, -14)

def tiny():
    """
    @rtype: float
    @return: 10^-20
    """
    return 1*pow(10, -20)

def timeFactor():
    """
    @rtype: float
    @return: Timestep scaling factor for propagation.
    """
    return 48.88821290839616

def invTimeFactor():
    """
    @rtype: float
    @return: Inverse of the timestep scaling factor for propagation.
    """
    return 0.02045482828087295

def periodicBoundaryTolerance():
    """
    @rtype: float
    @return: Size of buffer zone for the periodic cell.
    """
    return 3.0

def sqrtCoulombConstant():
    """
    @rtype: float
    @return: sqrt(k), for electrostatic energy (kq1q2/r).
    """
    return 18.2226123264

def pressureFactor():
    """
    @rtype: float
    
    """
    return 69478.0593635551

def pdbVelScalingFactor():
    """
    @rtype: float
    @return: Scaling factor for velocities from a PDB file.
    """
    return 20.45482706

def coulombFactor():
    """
    @rtype: float
    @return: Scaling factor to convert between units of Vm and C.
    """
    return pow(299792458.0, 2)*0.0000001

def electronCharge():
    """
    @rtype: float
    @return: Charge of an electron in Coulombs.
    """
    return 1.6021892*pow(10, -19)

def aaTom():
    """
    @rtype: float
    @return: Conversion from Angstroms to meter.
    """
    return 1*pow(10, 10)

def avogadro():
    """
    @rtype: float
    @return: Avogadro's number (atoms per mol).
    """
    return 6.022045*pow(10, 23)

def amuToKg():
    """
    @rtype: float
    @return: Scaling factor to convert AMU to SI kg.
    """
    return 1.6605655*pow(10, -27)

def kcalToJoule():
    """
    @rtype: float
    @return: Scaling factor to convert kcal to SI Joules.
    """
    return 1.0/4184.0

def fsTime():
    """
    @rtype: float
    @return: Number of fs per second (10^-15)
    """
    return 1*pow(10, 15)

def siBoltzmann():
    """
    @rtype: float
    @return: Boltzmann constant in SI units (J/K)
    """
    return 1.380662*pow(10, -23)
