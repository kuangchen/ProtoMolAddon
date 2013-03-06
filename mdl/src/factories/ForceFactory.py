import sys
import BondForce
import AngleForce
import DihedralForce
import HarmDihedralForce
import ImproperForce
import SimpleFullForce
import CutoffForce
import EwaldForce
import PMEForce 
import MultiGridForce   


class ForceFactory:
  def __init__(self):
    """
    Initializes mappings from boundary conditions (for bonded forces)
    and boundary conditions, algorithms and switching functions (for nonbonded
    forces) to SWIG-wrapped force object constructors (not instances, saving
    memory).
    """
    self.bondForces = {'Vacuum':BondForce.BSF_Vacuum,
                       'Periodic':BondForce.BSF_Periodic} #: Mapping from boundary conditions to two-atom bond force object constructors

    self.angleForces = {'Vacuum':AngleForce.ASF_Vacuum,
                        'Periodic':AngleForce.ASF_Periodic} #: Mapping from boundary conditions to three-atom angle force object constructors

    self.dihedralForces = {'Vacuum':DihedralForce.DSF_Vacuum,
                           'Periodic':DihedralForce.DSF_Periodic} #: Mapping from boundary conditions to four-atom dihedral force object constructors

    self.improperForces = {'Vacuum':ImproperForce.ISF_Vacuum,
                           'Periodic':ImproperForce.ISF_Periodic} #: Mapping from boundary conditions to four-atom improper force object constructors

    self.harmDihedralForces = {'Vacuum':HarmDihedralForce.HDSF_Vacuum,
                               'Periodic':HarmDihedralForce.HDSF_Periodic} #: Mapping from boundary conditions to harmonic dihedral force object constructors

    self.hd = 0  #: Number of harmonic dihedral forces
    
    ######################################################################################################
    # THE FOLLOWING NONBONDED FORCE OBJECTS ARE PREDEFINED AND WRAPPED AS SHARED OBJECTS FOR EFFICIENCY  #
    # OBJECTS ARE REGISTERED IN THE FACTORY AS MULTI-LEVEL MAPPINGS                                      #
    # INDEX ONE: BOUNDARY CONDITIONS                                                                     #
    # INDEX TWO: ALGORITHM                                                                               #
    # INDEX THREE: SWITCHING FUNCTION(S)                                                                 #                                  
    ######################################################################################################
    self.ljForces = {'Vacuum':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_V_U_L},
                               'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPVBC_C1SF_LJF,
                                         'C2':CutoffForce.NCSF_CCM_OAPVBC_C2SF_LJF,
                                         'Cutoff':CutoffForce.NCSF_CCM_OAPVBC_CSF_LJF,
                                         'Cn':CutoffForce.NCSF_CCM_OAPVBC_CNSF_LJF,
                                         'CmpCnCn':CutoffForce.NCSF_CCM_OAPVBC_CCNCNSF_LJF
                                         }
                               },
                     'Periodic':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_P_U_L},
                                 'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPPBC_C1SF_LJF,
                                           'C2':CutoffForce.NCSF_CCM_OAPPBC_C2SF_LJF,
                                           'Cutoff':CutoffForce.NCSF_CCM_OAPPBC_CSF_LJF,
                                           'Cn':CutoffForce.NCSF_CCM_OAPPBC_CNSF_LJF,
                                           'CmpCnCn':CutoffForce.NCSF_CCM_OAPPBC_CCNCNSF_LJF
                                           }
                                 }
                     }  #: Maps boundary conditions, algorithm and switching function to van der Waals force object constructor            

               
    self.cdeForces = {'Vacuum':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_V_U_CDE},
                                'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPVBC_C1SF_CFDE,
                                          'C2':CutoffForce.NCSF_CCM_OAPVBC_C2SF_CFDE,
                                          'Cutoff':CutoffForce.NCSF_CCM_OAPVBC_CSF_CFDE,
                                          'Cn':CutoffForce.NCSF_CCM_OAPVBC_CNSF_CFDE,
                                          'ComplementCn':CutoffForce.NCSF_CCM_OAPVBC_CMPCNNSF_CFDE}},
                      'Periodic':{'Cutoff':{'Cn':CutoffForce.NCSF_CCM_OAPPBC_CNSF_CFDE}}
                      } #: Maps boundary conditions, algorithm and switching function to coulomb dielectric force object constructor - used for implicit solvation
                                               

    self.coulombForces = {'Vacuum':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_V_U_C},
                                    'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPVBC_C1SF_CF,
                                              'C2':CutoffForce.NCSF_CCM_OAPVBC_C2SF_CF,
                                              'Cutoff':CutoffForce.NCSF_CCM_OAPVBC_CSF_CF,
                                              'Cn':CutoffForce.NCSF_CCM_OAPVBC_CNSF_CF,
                                              'ComplementCn':CutoffForce.NCSF_CCM_OAPVBC_CCNCNSF_CF},
                                    'SCPISM':{'C1':CutoffForce.NCSF_CCM_OAPVBC_C1SF_CSCPF,
                                              'C2':CutoffForce.NCSF_CCM_OAPVBC_C2SF_CSCPF,
                                              'Cn':CutoffForce.NCSF_CCM_OAPVBC_CNSF_CSCPF,
                                              'CmpCnCn':CutoffForce.NCSF_CCM_OAPVBC_CCNSF_CSCPF,
                                              'Cutoff':CutoffForce.NCSF_CCM_OAPVBC_CSF_CSCPF},
                                    'GB':{'Universal':CutoffForce.NCSF_CCM_OAPVBC_U_GB,
				          'C2':CutoffForce.NCSF_CCM_OAPVBC_C2_GB,
					  'Cn':CutoffForce.NCSF_CCM_OAPVBC_CN_GB},
			            'GBACE':{'Universal':CutoffForce.NCSF_CCM_OAPVBC_U_GBACE,
				             'C2':CutoffForce.NCSF_CCM_OAPVBC_C2_GBACE,
					     'Cn':CutoffForce.NCSF_CCM_OAPVBC_CN_GBACE},
				    'Ewald':{'Cutoff':EwaldForce.EWALD_V_TTT_CUTOFF,
                                             'C1':EwaldForce.EWALD_V_TTT_C1},
                                    'PME':{'Cutoff':PMEForce.PME_V_TTT_B,
                                           'C1':PMEForce.PME_V_TTT_B_C1},
                                    'MultiGrid':{'Cutoff':MultiGridForce.MG_V_TTT}},

                          'Periodic':{'SimpleFull':{'Universal':SimpleFullForce.NSFSF_P_U_C},
                                      'Cutoff':{'C1':CutoffForce.NCSF_CCM_OAPPBC_C1SF_CF,
                                                'C2':CutoffForce.NCSF_CCM_OAPPBC_C2SF_CF,
                                                'Cutoff':CutoffForce.NCSF_CCM_OAPPBC_CSF_CF,
                                                'Cn':CutoffForce.NCSF_CCM_OAPPBC_CNSF_CF,
                                                'ComplementCn':CutoffForce.NCSF_CCM_OAPPBC_CCNCNSF_CF},
                                      'Ewald':{'Cutoff':EwaldForce.EWALD_P_TTT_CUTOFF,
                                               'C1':EwaldForce.EWALD_P_TTT_C1},
                                      'PME':{'Cutoff':PMEForce.PME_P_TTT_B,
                                             'C1':PMEForce.PME_P_TTT_B_C1},
                                      'MultiGrid':{'Cutoff':MultiGridForce.MG_P_TTT}}
                                      

                          } #: Maps boundary conditions, algorithm and switching function to electrostatic force object constructor.  

    #self.bornForces = {'Vacuum':{'Cutoff':{'C1':CutoffForce.NCBF_CCM_OAPVBC_C1SF_CBF,
    #                                       'C2':CutoffForce.NCBF_CCM_OAPVBC_C2SF_CBF,
    #                                       'Cn':CutoffForce.NCBF_CCM_OAPVBC_CNSF_CBF,
    #                                       'CmpCnCn':CutoffForce.NCBF_CCM_OAPVBC_CCNSF_CBF,
    #                                       'Cutoff':CutoffForce.NCBF_CCM_OAPVBC_CSF_CBF}}} #: Maps boundary conditions, algorithm and switching function to an object which computes Born Radii. 

               
    self.ljCoulombForces = {'Vacuum':{'SimpleFull':{('Universal', 'Universal'):SimpleFullForce.NSFSF_V_U_L_U_C},
                                      'Cutoff':{('C2', 'C1'):CutoffForce.NCSF_CCM_OAPTVBC_C2SF_LJF_C1SF_CF,
                                                ('C2', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTVBC_C2SF_LJF_CSF_CF,
                                                ('Cn', 'Cn'):CutoffForce.NCSF_CCM_OAPTVBC_CN_LJF_CN_CF,
                                                ('Cn', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTVBC_CN_LJF_C_CF,
                                                ('CmpCnCn', 'CmpCnCn'):CutoffForce.NCSF_CCM_OAPTVBC_CN_LJF_CN_CF,
                                                ('CmpCnCn', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTVBC_CN_LJF_C_CF}

                                      },
                            'Periodic':{'SimpleFull':{('Universal', 'Universal'):SimpleFullForce.NSFSF_P_U_L_U_C},
                                        'Cutoff':{('C2', 'C1'):CutoffForce.NCSF_CCM_OAPTPBC_C2SF_LJF_C1SF_CF,
                                                  ('C2', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTPBC_C2SF_LJF_CSF_CF,
                                                  ('Cn', 'Cn'):CutoffForce.NCSF_CCM_OAPTPBC_CN_LJF_CN_CF,
                                                  ('Cn', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTPBC_CN_LJF_C_CF,
                                                  ('CmpCnCn', 'CmpCnCn'):CutoffForce.NCSF_CCM_OAPTPBC_CN_LJF_CN_CF,
                                                  ('CmpCnCn', 'Cutoff'):CutoffForce.NCSF_CCM_OAPTPBC_CN_LJF_C_CF}

                                        }
                            } #: Maps boundary conditions, algorithm and switching function pair to a unifed van der Waals and electrostatic force object constructor.  This saves performance by determining atom pairs just once for both types of pairwise forces.


               
  def createBondForce(self,bc):
    """
    Return a bond force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped bond force object.
    """
    return self.bondForces[bc]()
      
  def createAngleForce(self, bc):
    """
    Return an angle force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped angle force object.
    """
    return self.angleForces[bc]()
   
  def createDihedralForce(self, bc):
    """
    Return a dihedral force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped dihedral force object.
    """
    return self.dihedralForces[bc]()

  def createImproperForce(self, bc):
    """
    Return an improper force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @rtype: Force
    @return: SWIG-wrapped improper force object.
    """
    return self.improperForces[bc]()

  def createHarmDihedralForce(self, bc, params):
    """
    Return a harmonic dihedral force object.

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @type params: dict
    @param params: Mapping from parameter names to corresponding values.

    @rtype: Force
    @return: SWIG-wrapped harmonic dihedral force object.
    """
    if (str(type(params['kbias']))[7:11] == 'list'):
      self.hd += 1  # INCREASE HARMONIC DIHEDRAL COUNTER
      return self.harmDihedralForces[bc](params['kbias'][self.hd-1], params['dihedralnum'][self.hd-1], params['angle'][self.hd-1], 0)
    else:
      return self.harmDihedralForces[bc](params['kbias'], params['dihedralnum'], params['angle'], 0)      

  def createLennardJonesForce(self, bc, params):
      """
      Return a van der Waals force object.
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)
      
      @type params: dict
      @param params: Mapping from parameter names to corresponding values.

      @rtype: Force
      @return: SWIG-wrapped van der Waals force object.
      """
      alg = self.getParameter(params, 'algorithm', "SimpleFull")
      switch = self.getParameter(params, 'switching', "Universal")
      newforce = self.lookup(self.ljForces, 'LennardJones', bc, alg, switch)
      print newforce
      return self.applyParameters(newforce, bc, alg, switch, params)

  def createCoulombDiElecForce(self, bc, params):
      """
      Return an coulomb dielectric force object (for implicit solvation).
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)
      
      @type params: dict
      @param params: Mapping from parameter names to corresponding values.

      @rtype: Force
      @return: SWIG-wrapped coulomb dielectric force object.
      """
      alg = self.getParameter(params, 'algorithm', "SimpleFull")
      switch = self.getParameter(params, 'switching', "Universal")
      newforce = self.lookup(self.cdeForces, 'CoulombDiElec', bc, alg, switch)
      return self.applyParameters(newforce, bc, alg, switch, params)
       	
  def createCoulombForce(self, bc, params):
      """
      Return an electrostatic force object.
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)
      
      @type params: dict
      @param params: Mapping from parameter names to corresponding values.

      @rtype: Force
      @return: SWIG-wrapped electrostatic force object.
      """
      alg = self.getParameter(params, 'algorithm', "SimpleFull")
      switch = self.getParameter(params, 'switching', "Universal")
      newforce = self.lookup(self.coulombForces, 'Coulomb', bc, alg, switch)
      return self.applyParameters(newforce, bc, alg, switch, params)
      
  def createBornBurial(self):
      #return SimpleFullForce.NSFSF_V_U_GB_GBORNBUR()
      return SimpleFullForce.NSFSF_V_U_GB_GBORNBUR().makeNew()

  def createBornRadii(self):
      #return SimpleFullForce.NSFSF_V_U_GB_GBORN()
      return SimpleFullForce.NSFSF_V_U_GB_GBORN().makeNew()
      
  def createBornForce(self, bc, params):
      """
      Return an electrostatic force object.
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)
      
      @type params: dict
      @param params: Mapping from parameter names to corresponding values.

      @rtype: Force
      @return: SWIG-wrapped force object which computes Born Radii.
      """
      # Modified TMC 6/30/08: ALWAYS use Cutoff switch and alg
      alg = "Cutoff"
      newforce = self.bornForces[bc][alg]['Cutoff']()
      return self.applyParameters(newforce, bc, alg, 'Cutoff', params)

  def createLennardJonesCoulombForce(self, bc, params):
      """
      Return a coupled van der Waals - electrostatic force evaluation object.
      Saves performance by calculating atom pairs only once for both
      type of forces.
      
      @type bc: string
      @param bc: Boundary conditions (Periodic or Vacuum)
      
      @type params: dict
      @param params: Mapping from parameter names to corresponding values.

      @rtype: Force
      @return: SWIG-wrapped coupled van der Waals and electrostatic force object.
      """
      alg = self.getParameter(params, 'algorithm', "Cutoff")  # Different default
      switch = self.switchingFunctions(params, 2)
      newforce = self.ljCoulombForces[bc][alg][switch]()
      return self.applyParameters(newforce, bc, alg, switch, params)

  def switchingFunctions(self, params, number):
      """
      Little helper function, to take a set of parameters
      and obtain a specific number of switching functions,
      setting defaults as necessary (i.e. if no switching functions
      are provided, Universal is assumed)

      @type params: dict
      @param params: Mapping from parameter names to values

      @type number: int
      @param number: Number of switching functions to obtain
      """
      fxns = ()
      if (not params.has_key('switching')):
        for i in range(0, number):
           fxns += ("Universal",)
      else:
        if (number == 1):
           fxns += (params['switching'],)
        else:
           if (str(type(params['switching']))[7:11] == 'list'):
             for j in range(0, len(params['switching'])):
                fxns += (params['switching'][j],)
             for k in range(len(params['switching']), number):
                fxns += (fxns[k-1],)
           else:
             for j in range(0, number):
               fxns += (params['switching'],)  # Just add the first
      return fxns


  def lookup(self, table, forcename, bc, alg, switch):
    """
    Helper function, looks up the passed boundary conditions, algorithm
    and switching function in the passed table.  If it exists, return the
    corresponding force object instance.  Otherwise, output a helpful error
    message for the user.

    @type table: dict
    @param table: Mapping from parameters to force constructor.

    @type forcename: string
    @param forcename: Name of force
    
    @type bc: string
    @param bc: Boundary conditions

    @type alg: string
    @param alg: Nonbonded calculation algorithm.

    @type switch: string or tuple of strings
    @param switch: Switching function(s).
    """

    if (not table.has_key(bc)):
      print "[MDL] ERROR: UNRECOGNIZED BOUNDARY CONDITIONS: ", bc, " FOR FORCE ", forcename
      print "POSSIBILITIES INCLUDE: ",
      for i in table.iterkeys():
        print i
    else:
      if (not table[bc].has_key(alg)):
        print "[MDL] ERROR: UNRECOGNIZED ALGORITHM: ", alg,
        print " FOR FORCE ", forcename, " WITH ", bc, " BOUNDARY CONDITIONS."
        print "POSSIBILITIES INCLUDE: ",
        for i in table[bc].iterkeys():
          print i
      else:
        if (not table[bc][alg].has_key(switch)):
          print "[MDL] ERROR: UNRECOGNIZED SWITCHING FUNCTION(S): ", switch,
          print " FOR FORCE ", forcename, " WITH ", bc, "BOUNDARY CONDITIONS",
          print " AND ", alg, " ALGORITHM."
          print "POSSIBILITIES INCLUDE: ",
          for i in table[bc].iterkeys():
            print i
        else:
          return table[bc][alg][switch]()

    
    
    
  # Return the parameter if it exists, otherwise return the provided default value
  def getParameter(self, params, name, defaultval=0):
    """
    Return a parameter value if it exists, otherwise return the passed default value.

    @type params: dict
    @param params: Mapping from parameter names to values

    @type name: string
    @param name: Name of the parameter

    @type defaultval: (any type)
    @param defaultval: Default value for the parameter if it is not found
    """
    if (params.has_key(name)):
      return params[name]
    else:
      return defaultval
    
  # Add application of parameters here
  # This is its own function, because a lot of algorithm-switching function combinations
  # will have the same extra parameters, and this avoids code repetition
  def applyParameters(self, newforce, bc, alg, switch, params, fastelectro=None):
    """
    Depending on the algorithm used, set parameter values for the passed
    force object.

    @type newforce: Force
    @param newforce: Pairwise force object

    @type bc: string
    @param bc: Boundary conditions (Periodic or Vacuum)

    @type alg: string
    @param alg: Algorithm for pairwise evaluation

    @type switch: string
    @param switch: Switching function

    @type params: dict
    @param params: Mapping from parameter names to values

    @type fastelectro: dict
    @param fastelectro: Mapping from fast electrostatic parameter names to values

    @rtype: Force
    @return: Pairwise force object with instantiated parameters.
    """
    if (str(newforce).find('OneAtomPairTwo') != -1): # WE HAVE A LJ_COU PAIR
      if (alg == "Cutoff"):
        newforce = newforce.makeNewPair(switch[0], switch[1],
                                        self.getParameter(params,'cutoff'),
                                        self.getParameter(params,'switchon',0),
                                        self.getParameter(params,'switchoff',self.getParameter(params,'cutoff')),
                                        self.getParameter(params,'order',2))
      else:
        newforce = newforce.makeNewPair(self.getParameter(params, 'blocksize', 32))
    else:
      if (str(newforce).find('CoulombForceDiElec') != -1): # Dielectric forces have extra params
        eps = self.getParameter(params, 'epsilon', 1)
        d = self.getParameter(params, 'd', 78)
        s = self.getParameter(params, 's', 0.3)
        if (alg == "SimpleFull"):
          newforce = newforce.makeNewDiElec(self.getParameter(params, 'blocksize', 32), eps, d, s)
        elif (alg == "Cutoff"):
          if (switch == "C1"):
            newforce = newforce.makeNewDiElec(self.getParameter(params, 'cutoff'), eps, d, s)
          elif (switch == "C2"):
            newforce = newforce.makeNewDiElec(self.getParameter(params, 'cutoff'),
                                              eps, d, s,
                                              self.getParameter(params, 'switchon'))
          else:
            newforce = newforce.makeNewDiElec(self.getParameter(params, 'cutoff'),
                                              eps, d, s,
                                              self.getParameter(params, 'switchon'),
                                              self.getParameter(params, 'switchoff'),
                                              self.getParameter(params, 'order', 2))      
      elif (alg == "Ewald"):
         alpha = self.getParameter(params, 'alpha', -1)
         accuracy = self.getParameter(params, 'accuracy', 0)
         expansionfactor = self.getParameter(params, 'expansionFactor', 0)
         if (bc == "Periodic"):
            newforce = newforce.makeNew(alpha, accuracy)
         else:
            newforce = newforce.makeNew(alpha, accuracy, expansionfactor)
      elif (alg == "PME"):
         alpha = self.getParameter(params, 'alpha', -1)
         accuracy = self.getParameter(params, 'accuracy', 0.000001)
         order = self.getParameter(params, 'order', 4)
         cutoff = self.getParameter(params, 'cutoff')
         gridsize = self.getParameter(params, 'gridsize')
         expfactor = self.getParameter(params, 'expfactor', 3)
         if (bc == "Periodic"):
            newforce = newforce.makeNew(gridsize, cutoff, order, alpha, accuracy)
         else:
            newforce = newforce.makeNew(gridsize, cutoff, order, alpha, accuracy, expfactor)
      elif (alg == "MultiGrid"):
         s = self.getParameter(params, 'smoothdist')
         if (bc == "Periodic"):
            gridsize = self.getParameter(params, 'gridsize')
         else:
            h = self.getParameter(params, 'h', 4)
            o = self.getParameter(params, 'o', 0)
         levels = self.getParameter(params, 'levels')
         order = self.getParameter(params, 'order', 4)
         ratio = self.getParameter(params, 'ratio', 2)
         if (bc == "Periodic"):
            newforce = newforce.makeNew(s, gridsize, levels, order, ratio)
         else:
            newforce = newforce.makeNew(s, h, o, levels, order, ratio)
      else:
        if (alg == "SimpleFull"):
          newforce = newforce.makeNew(self.getParameter(params, 'blocksize', 32))
        elif (alg == "Cutoff"):
          if (str(newforce).find('Born') != -1):
            newforce = newforce.makeNew(self.getParameter(params, 'borncutoff'))
          else:
            if (switch == "C1" or switch == "Cutoff"):
               newforce = newforce.makeNew(self.getParameter(params, 'cutoff'))
            elif (switch == "C2"):
              newforce = newforce.makeNew(self.getParameter(params, 'cutoff'),
                                          self.getParameter(params, 'switchon'))
            else:
              newforce = newforce.makeNew(self.getParameter(params, 'cutoff'),
                                          self.getParameter(params, 'switchon'),
                                          self.getParameter(params, 'switchoff'),
                                          self.getParameter(params, 'order', 2))
        elif (alg == "GB"):
	   if (switch == "Universal"):
	      newforce = newforce.makeNewGB(self.getParameter(params, 'solute', 1.0),
	                                    self.getParameter(params, 'solvent', 80.0))
	   elif (switch == "C2"):
	      newforce = newforce.makeNewGB(self.getParameter(params, 'solute', 1.0),
	                                    self.getParameter(params, 'solvent', 80.0),
					    self.getParameter(params, 'cutoff'),
					    self.getParameter(params, 'switchon'))
           else:
              newforce = newforce.makeNewGB(self.getParameter(params, 'solute', 1.0),
	                                    self.getParameter(params, 'solvent', 80.0),
					    self.getParameter(params, 'cutoff'),
					    self.getParameter(params, 'switchon'),
					    self.getParameter(params, 'switchoff'),
					    self.getParameter(params, 'order', 2))
        elif (alg == "GBACE"):
	   newforce = SimpleFullForce.NSFSF_V_U_GB_GB().makeNew()
	   #if (switch == "Universal"):
	   #   newforce = newforce.makeNewGB(self.getParameter(params, 'solvation', 2.26/418.4),
	   #                                 self.getParameter(params, 'sphere', 1.4))
	   #elif (switch == "C2"):
	   #   newforce = newforce.makeNewGB(self.getParameter(params, 'solvation', 2.26/418.4),
	   #                                 self.getParameter(params, 'sphere', 1.4),
	   #				    self.getParameter(params, 'cutoff'),
	   #				    self.getParameter(params, 'switchon'))
           #else:
           #   newforce = newforce.makeNewGB(self.getParameter(params, 'solvation', 2.26/418.4),
	   #                                 self.getParameter(params, 'sphere', 1.4),
	   #				    self.getParameter(params, 'cutoff'),
	   #				    self.getParameter(params, 'switchon'),
	   #				    self.getParameter(params, 'switchoff'),
	   #				    self.getParameter(params, 'order', 2))
 
    return newforce
 
    return newforce
