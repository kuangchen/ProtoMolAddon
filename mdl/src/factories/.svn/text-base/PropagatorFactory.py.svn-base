import sys
import types
import os
import Constants
import TopologyUtilities
import copy
import numpy
#import dl
#import os
#if os.name == 'posix':
#   sys.setdlopenflags(dl.RTLD_LAZY|dl.RTLD_GLOBAL)


# THIS SUBROUTINE INITIALIZES, SETS UP AND RUNS PRE-
# AND POST- INITIALIZATION MODIFIERS ON THE PASSED PROPAGATOR.
# ALSO UPDATES SYSTEM CENTER OF MASS AND MOMENTA.
# USER SHOULD NOT NEED TO INVOKE THIS DIRECTLY.
def setPropagator(prop, phys, forces, obj, levelswitch=False):
	"""
	Set and initialize a propagator object.

	@type prop: Propagator
	@param prop: MDL Propagator object

	@type phys: Physical
	@param phys: The physical system.

	@type forces: Forces
	@param forces: MDL Forces object
	
	@type obj: STS/MTS
	@param obj: Prototyped propagator object

	@type levelswitch: boolean
	@param levelswitch: True if we are changing levels in the hierarchy.  Default false.
	"""
        if (prop.myLevel == 0):
           forces.forcevec = obj.getForces()
	   if (dir(obj).count('setIntegratorSetPointers') != 0):
		   obj.setIntegratorSetPointers(obj, phys.myEig, 1)
           if (not hasattr(phys, 'app')):
	      phys.app = obj.appInit(phys.myTop,phys.posvec,phys.velvec,forces.energies)
	   else:
	      obj.initialize(phys.app)
           phys.app.energies = forces.energies
           #obj.initialize(phys.myTop,phys.posvec,phys.velvec,forces.energies)
        # Do not perform garbage collection if we are setting our propagator
	# to something else, and aren't simply changing levels in the hierarchy
	
	if (prop.myPropagator != 0 and (not levelswitch)):
	    del(prop.myPropagator)
	
	prop.myPropagator = obj
	if (prop.isMDL(prop.myPropagator)):
           prop.runModifiers(prop.myPropagator.preinitmodifiers, phys, forces, prop, prop.myPropagator)
           prop.myPropagator.init(phys, forces, prop)
           prop.runModifiers(prop.myPropagator.postinitmodifiers, phys, forces, prop, prop.myPropagator)
        phys.updateCOM_Momenta()

# RUN THE PASSED PROPAGATOR FOR numsteps
def executePropagator(prop, phys, forces, io, numsteps):
      """
      Run and finish the propagator.
      
      @type prop: Propagator
      @param prop: MDL Propagator object
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object

      @type io: IO
      @param io: MDL IO object
      
      @type numsteps: int
      @param numsteps: Number of steps to run
      """
      # RUN
      if (prop.myStep == 0 and prop.myLevel == 0):
           io.run(phys, forces, prop.myStep, prop.myTimestep, prop.myPropagator)
      if (prop.isMDL(prop.myPropagator)):
           prop.runModifiers(prop.myPropagator.prerunmodifiers, phys, forces, prop, prop.myPropagator)
      ii = 0

      #prop.myPropagator.run(numsteps)
      #return
      while (ii < numsteps):
              nextstop = int(numpy.minimum(numsteps, io.computeNext(ii, phys.remcom, phys.remang)))

              if (prop.isMDL(prop.myPropagator)):
                 while (ii < nextstop):
		    prop.myPropagator.run(phys, forces, prop)
		    ii += 1
                 phys.myTop.time = prop.myStep*prop.myPropagator.getTimestep()*Constants.invTimeFactor()
                 phys.updateCOM_Momenta()
              else:
                 prop.myPropagator.run(nextstop-ii)
                 ii = nextstop
              
	      
	      if (prop.myLevel == 0):
                 prop.myStep = nextstop
		 phys.myTop.time = prop.myStep*prop.myPropagator.getTimestep()
		 io.run(phys, forces, prop.myStep, prop.myTimestep, prop.myPropagator)
	      if (phys.remcom > 0 and ii % phys.remcom == 0):
                 TopologyUtilities.remLin(phys.velvec, phys.myTop)
              if (phys.remang > 0 and ii % phys.remang >= 0):
	         TopologyUtilities.remAng(phys.posvec, phys.velvec, phys.myTop)
      if (prop.isMDL(prop.myPropagator)):
           prop.runModifiers(prop.myPropagator.postrunmodifiers, phys, forces, prop, prop.myPropagator)
           prop.myPropagator.finish(phys, forces, prop)




def _get_mod(modulePath):
    """
    Given a path to a Python module, return the module itself.
    This allows the factory to search for example the propagators/ and
    modifiers/ packages

    @type modulePath: string
    @param modulePath: Path to the Python module.

    @rtype: Python module
    @return: The Python module.
    """
    try:
        aMod = sys.modules[modulePath]
        if not isinstance(aMod, types.ModuleType):
            raise KeyError
    except KeyError:
        # The last [''] is very important!
        aMod = __import__(modulePath, globals(), locals(), [''])
        sys.modules[modulePath] = aMod
    return aMod


class PropagatorFactory:
   """
   Contains mappings from propagator names to creation functions,
   their types (object or method) and parameters with default values.
   Given a name and parameter values, the propagator factory can
   return either a method handle or an instance of a propagator object.
   """
   def __init__(self):
      #self.objects = []  #: A list of propagator objects created.
      
      self.registry = dict()  #: Maps propagator names to creation functions and parameters/defaults
      self.modifier_registry = dict()  #: Maps modifier names to method handles

      self.registerAllPM()
      # REGISTER ALL MDL PROPAGATORS HERE.
      # TO DO THIS, WE MUST FIRST OBTAIN
      # A LIST OF ALL FILES IN propagators/objects:
      mdlpath = os.getenv('MDLROOT')
      myList = os.listdir(mdlpath+'/src/propagators/classes')

      # NOW REMOVE THE .svn/ directory
      if (myList.count('.svn') != 0): myList.remove('.svn')

      # NOW PROCESS THE LIST.
      # 1. REMOVE DIRECTORIES
      # 2. REMOVE NON .py FILES
      i = 0
      while (i < myList.__len__()):
          if ( (not os.path.isfile(mdlpath+'/src/propagators/classes/'+myList[i])) or
               (not myList[i].endswith('.py')) ):
              myList.remove(myList[i])
          else:
              i += 1

      # NOW WE ARE LEFT WITH .py FILES.
      # WE MUST STRIP THIS ENDING TO GET MODULE NAMES
      for i in range(0, myList.__len__()):
          myList[i] = myList[i][0:myList[i].__len__()-3]

      # NOW, REMOVE ALL MODULES THAT ARE ALREADY REGISTERED
      i = 0
      while i < len(myList):
          if (self.registry.has_key(myList[i])):
              print "[MDL] WARNING: REMOVING SECOND ", myList[i]
              myList.remove(myList[i])
          else:
              i += 1

      # MDL PROPAGATORS ARE LEFT
      ## NOW REGISTER.
      for i in range(0, myList.__len__()):
          mod = _get_mod('classes.'+myList[i])
          if (hasattr(mod, 'name') and hasattr(mod, 'parameters')):
              self.registerObject(mod, mod.name, mod.parameters)
          if (hasattr(mod, 'modifiers')):
              for j in range(0, mod.modifiers.__len__()):
                  self.registerModifier(mod.modifiers[j][0],
                                        mod.modifiers[j][1],
                                        mod.name)

      # REPEAT THE PROCESS FOR THE methods DIRECTORY
      # EXCEPT REGISTER METHODS.
      myList = os.listdir(mdlpath+'/src/propagators/functions')
      if (myList.count('.svn') != 0): myList.remove('.svn')

      i = 0
      while (i < myList.__len__()):
          if ( (not os.path.isfile(mdlpath+'/src/propagators/functions/'+myList[i])) or
               (not myList[i].endswith('.py')) ):
              myList.remove(myList[i])
          else:
              i += 1
              
      for i in range(0, myList.__len__()):
          myList[i] = myList[i][0:myList[i].__len__()-3]
      i = 0
      while i < len(myList):
          if (self.registry.has_key(myList[i])):
              print "[MDL] WARNING: REMOVING SECOND ", myList[i]
              myList.remove(myList[i])
          else:
              i += 1
              
      for i in range(0, myList.__len__()):
          mod = _get_mod('functions.'+myList[i])
          if (hasattr(mod, 'name') and hasattr(mod, 'parameters')):
              self.registerMethod(mod, mod.name, mod.parameters)
          if (hasattr(mod, 'modifiers')):
              for j in range(0, mod.modifiers.__len__()):
                  self.registerModifier(mod.modifiers[j][0],
                                        mod.modifiers[j][1],
                                        name)


   def query(self, name="None"):
      """
      Query the factory for information on all available propagation schemes
      This will display their type (object or method) and parameters and
      default values.
      You can also query with a name to get information about a specific
      propagation scheme

      @type name: string
      @param name: Propagation scheme name, default is None
      """
      print "***********************************************************"
      print "*             MDL AVAILABLE PROPAGATORS                   *"
      print "***********************************************************"
      print ""
      print "Name".ljust(25),
      print "Type".ljust(10),
      print "Parameters/Defaults"
      print ""
      if (name != "None"):
              if (not self.registry.has_key(name)):
                 print "Propagator ", name, " is not defined."
              else: 
                 print name.ljust(25),
                 print self.registry[name]['type'].ljust(10),
                 j = 0
                 while (j < len(self.registry[name]['defaults'])):
                    print self.registry[name]['defaults'][j],
                    print ':',
                    print self.registry[name]['defaults'][j+1]
                    print "".ljust(25),
                    print "".ljust(10),
                    j += 2
      else:
              for i in self.registry.keys():
                 print str(i).ljust(25),
                 print self.registry[i]['type'].ljust(10),
                 j = 0
                 while (j < len(self.registry[i]['defaults'])):
                    print self.registry[i]['defaults'][j],
                    print ':',
                    print self.registry[i]['defaults'][j+1]
                    print "".ljust(25),
                    print "".ljust(10),
                    j += 2
                 print ""

   # Name: Integrator name, string
   # Defaults: Parameters and values, tuple
   def registerObject(self, myModule, name, defaults):
      """
      Register a propagator object prototyped in Python

      @type myModule: Python Module
      @param myModule: Module in which the class definition is contained

      @type name: string
      @param name: Name of propagation scheme

      @type defaults: tuple
      @param defaults: Mapping from parameter names to default values
      """
      self.registry[name] = {'constructor':getattr(myModule, name),
                             'type':'object',
                             'defaults':defaults}

   def registerMethod(self, myModule, name, defaults):
      """
      Register a propagator Python method.

      @type myModule: Python Module
      @param myModule: Module in which the method definition is contained

      @type name: string
      @param name: Name of propagation scheme

      @type defaults: tuple
      @param defaults: Mapping from parameter names to default values
      """
      self.registry[name] = {'constructor':getattr(myModule, name),
                             'type':'method',
                             'defaults':defaults}

   def registerPMObject(self, name, defaults):
      """
      Register a SWIG-wrapped propagator object

      @type name: string
      @param name: Name of propagation scheme

      @type defaults: tuple
      @param defaults: Mapping from parameter names to default values
      """
      # This is what happens when naming conventions are not followed....
      if (name[0:4] != "Norm" and name[0:4] != "Nume" and name != "LeapfrogTruncatedShadow" and name != "HessianInt"):
	      modname = name+"Integrator"
      else:
	      modname = name
      myModule = _get_mod(modname)   # For all SWIG-wrapped, the module name is the same as the propagator name
      self.registry[name] = {'constructor':getattr(myModule, modname),
                             'type':'protomol',
                             'defaults':defaults}


   def registerModifier(self, name, thetype, theprop):
      """
      Register a propagator modifier

      @type name: string
      @param name: Name of propagation scheme

      @type thetype: string
      @param thetype: Type of modifier (preinit, postinit, preforce, ...)
      
      @type theprop: string
      @param theprop: Name of the propagator object that this modifier operates upon.
      """
      myModule = _get_mod('modifiers.'+name)
      self.modifier_registry[name] = {'constructor':getattr(myModule, name),
                                      'type':thetype,
                                      'prop':theprop}

   def getType(self, name):
      """
      Return the type (class or function) of a propagator

      @type name: string
      @param name: Name of propagation scheme

      @rtype: string
      @return: Type of the propagator
      """
      if (not self.registry.has_key(name)):
	      print "[MDL] ERROR: NO PROPAGATOR NAMED", name
              sys.exit(1)
      return self.registry[name]['type']

   def findArg(self, name, pars):
      """
      Find a parameter name in a parameter tuple.

      @type name: string
      @param name: Name of parameter

      @type pars: tuple
      @param pars: Parameter list

      @rtype: int
      @return: Index of the parameter in the list (-1 if not found)
      """
      ii = 0
      while ((ii < pars.__len__()) and (type(pars[ii]).__name__ == 'tuple')):
          if (pars[ii][0] == name):
              return ii
          ii += 1
      return -1

   def applyModifiers(self, obj, name):
       """
       Instantiate modifiers for a propagator.

       @type obj: STS/MTS
       @param obj: Propagator object

       @type name: string
       @param name: Propagator name

       @rtype: STS/MTS
       @return: Propagator object with instantiated modifiers.
       """
       for modname in self.modifier_registry.keys():
           entry = self.modifier_registry[modname]
           if (entry['prop'] == name):
                 if (entry['type'] == "PreInit"):
                     obj.preinitmodifiers.append(entry['constructor'])
                 if (entry['type'] == "PostInit"):
                     obj.postinitmodifiers.append(entry['constructor'])
                 if (entry['type'] == "PreForce"):
                     obj.preforcemodifiers.append(entry['constructor'])
                 if (entry['type'] == "PostForce"):
                     obj.postforcemodifiers.append(entry['constructor'])
                 if (entry['type'] == "PreRun"):
                     obj.prerunmodifiers.append(entry['constructor'])
                 if (entry['type'] == "PostRun"):
                     obj.postrunmodifiers.append(entry['constructor'])
       return obj
           
   def create(self, level, *args):
      """
      Accept a Python tuple containing the propagator name,
      timestep, parameter values and force field.
      This tuple can thus be various sizes depending on the number
      of parameters.
      Create and return a corresponding instantiated propagator object,
      or method handle.

      @type args: tuple
      @param args: List of propagator name, dt, parameter values and force field.  If the propagation scheme is MTS, this list is followed by the name of the next propagator in the chain and its associated parameters, etc.

      @rtype: STS, MTS, or Python method handle
      @return: The propagator.
      """
      global objects 
      if (level == 1): objects = []
      if (not self.registry.has_key(args[0])):
         if (type(args[0]).__name__ == 'str'):
             print "[MDL] ERROR: UNRECOGNIZED PROPAGATOR: ",
             print args[0]
             return
         if (args.__len__() <= 2):
            print "[MDL] WARNING: USING LEAPFROG AS INNER PROPAGATOR"
            return self.create(level, "Leapfrog", args[0], args[1])
         else:
            print "[MDL] WARNING: USING LEAPFROG AS INNER PROPAGATOR"
            return self.create(level, "Leapfrog", args[0], args[1], args.__getslice__(2, args.__len__()))
      regprop = self.registry[args[0]]
      if (regprop['type'] == "protomol"):
         arglist = (args[1],)  # Timestep or cyclelength
         ii = 0
         while (ii < regprop['defaults'].__len__()):
            if (args[3].has_key(regprop['defaults'][ii])):
                arglist += (args[3][regprop['defaults'][ii]],)
            else:
                arglist += (regprop['defaults'][ii+1],)
            ii += 2                     
	 arglist += (args[2],)
         if (args.__len__() > 4):
	    # Python referencing - we must save the object we create
	    # in a temporary array; otherwise when it loses scope it gets
	    # destroyed
            objects.append(apply(self.create, (level+1,)+args.__getslice__(4, len(args))))
	    objects[len(objects)-1].thisown=0
	    # Now append to the list
	    arglist += (objects[len(objects)-1],)
	 if (args[0] == "OpenMM"):
	    retval = apply(regprop['constructor'], ())
	    return retval.setParams(*arglist)
	 else:
            return apply(regprop['constructor'], arglist)
      elif (regprop['type'] == "object"):
         if (len(args) <= 4):
            objects.append(apply(regprop['constructor'], (args[1], args[2])))
	    objects[len(objects)-1].dt = args[1]*Constants.invTimeFactor()
         else:
            objects.append(apply(self.create, (level+1,)+args.__getslice__(4, args.__len__())))
            objects.append(apply(regprop['constructor'], (args[1], args[2], objects[objects.__len__()-1])))
	    objects[len(objects)-1].thisown=0
	    objects[len(objects)-1].cyclelength = args[1]
            objects[objects.__len__()-1].next = objects[objects.__len__()-2]
            if (hasattr(objects[objects.__len__()-1], "setMOLLYForceGroups")):
                objects[objects.__len__()-1].setMOLLYForceGroups()
         ii = 0
         while (ii < regprop['defaults'].__len__()):
            if (args[3].has_key(regprop['defaults'][ii])):
                setattr(objects[objects.__len__()-1], regprop['defaults'][ii], args[3][regprop['defaults'][ii]])
            else:
                setattr(objects[objects.__len__()-1], regprop['defaults'][ii], regprop['defaults'][ii+1])
            ii += 2
         objects[objects.__len__()-1].preinitmodifiers = list()
         objects[objects.__len__()-1].postinitmodifiers = list()
         objects[objects.__len__()-1].prerunmodifiers = list()
         objects[objects.__len__()-1].preforcemodifiers = list()
         objects[objects.__len__()-1].postforcemodifiers = list()
         objects[objects.__len__()-1].postrunmodifiers = list()
         return objects[objects.__len__()-1]
      elif (regprop['type'] == "method"):
         tup = (args[1], args[2], args[3], args[4], args[5], args[6])
         ii = 0
         while (ii < regprop['defaults'].__len__()):
             if (args[7].has_key(regprop['defaults'][ii])):
                 tup += (args[7][regprop['defaults'][ii]],)
             else:
                 tup += (regprop['defaults'][ii+1],)
             ii += 2

         ourtup = args
         i = 8;
         # LOOP WHILE THERE ARE MORE METHODS IN THE CHAIN
         while (len(ourtup) > i):
             tup += (self.registry[ourtup[i]]['constructor'],)  # NEXT METHOD
             ii = 0
             while (ii < self.registry[ourtup[i]]['defaults'].__len__()):
                 if (ourtup[i+3].has_key(self.registry[ourtup[8]]['defaults'][ii])):
                     tup += (ourtup[i+3][ii],)
                 else:
                     tup += (self.registry[ourtup[i]]['defaults'][ii+1],)
                 ii += 2
             tup += ourtup.__getslice__(i+1,i+3) # dt, ff
             ourtup = ourtup.__getslice__(i+4, len(ourtup))
             i = 4
         objects.append(apply(regprop['constructor'], tup))
         return objects[objects.__len__()-1]
                          
   def registerAllPM(self):
      """
      Register all SWIG-wrapped propagators.
      """
      self.registerPMObject("BBK",
                       ('temp', 300,
                        'gamma', 2,
                        'seed', 1234))
      self.registerPMObject("Impulse", ())
      self.registerPMObject("LangevinImpulse",
                       ('temp', 300,
                        'gamma', 2,
                        'seed', 1234))
      self.registerPMObject("LangevinLeapfrog",
                       ('temp', 300,
                        'gamma', 2,
                        'seed', 1234))
      self.registerPMObject("CGMinimizer",
      		       ('alpha', 0.001,
			'beta', 0.05,
			'restart', 0))
      self.registerPMObject("Leapfrog", ())
      self.registerPMObject("LeapfrogTruncatedShadow", ())
      self.registerPMObject("HessianInt",
			    ('eigvecFile', '',
			     'eigvalFile', '',
			     'hessianFile', '',
			     'sortByAbs', 0,
			     'fixedModes', 0,
			     'textEigFile', 0,
			     'fullDiag', 0,
                             'eigenValueThresh', 5.0,
                             'blockVectorCols', 0,
                             'residuesPerBlock', 1,
                             'blockCutoffDistance', 10,
                             'massweight', 1,
                             'noseMass', 0,
                             'autoParameters', 0,
                             'geometric', 0,
                             'numerically', 0))
      self.registerPMObject("NormalModeBrownian",
			    ('firstmode', -1,
			     'numbermodes', -1,
			     'gamma', 80,
			     'seed', 1234,
			     'temperature', 300.0,
			     'avForceFile', "",
			     'inForceFile', ""))
      self.registerPMObject("NormalModeDiagonalize",
			    ('reDiagFrequency', 0,
                             'fullDiag', 1,
                             'removeRand', 0,
                             'rediagHysteresis', 0,
                             'eigenValueThresh', 5,
                             'blockVectorCols', 0,
                             'residuesPerBlock', 1,
                             'blockCutoffDistance', 10,
			     'autoParameters', True,
			     'autoDT', False,
			     'pdm', False,
			     'minimLim', 0.1,
                             'maxminsteps', 100))
                            #('averageSteps', 1,
			    # 'avStepSize', 1.0,
			    # 'reDiagFrequency', 0,
			    # 'fullDiag', 0,
			    # 'removeRand', 0,
			    # 'minSteps', 0,
			    # 'minLim', 1.0))
      self.registerPMObject("NormalModeLangevin",
			    ('firstmode', 1,
			     'numbermodes', 1,
			     'gamma', 80,
			     'seed', 1234,
			     'temperature', 300,
			     'gencompnoise', 0))
      self.registerPMObject("NormalModeLangLf",
			    ('firstmode', 1,
			     'numbermodes', 1,
			     'gamma', 80,
			     'seed', 1234,
			     'temperature', 300,
			     'gencompnoise', 0))
      self.registerPMObject("NormalModeMinimizer",
			    ('firstmode', -1,
			     'numbermodes', -1,
			     'gamma', 80,
			     'seed', 1234,
			     'temperature', 300.0,
			     'minimlim', 0.1,
			     'randforce', 1,
			     'rediag', 0,
			     'simplemin', 1,
			     'redmaxmin', 0))
      self.registerPMObject("NormalModeMori",
			    ('firstmode', 1,
			     'numbermodes', 1,
			     'gamma', 80,
			     'seed', 1234,
			     'temperature', 300,
			     'modeOutput', "",
			     'instForce', 0))
      self.registerPMObject("NormalModeRelax",
			    ('minimlim', 0.1,
			     'rediag', 0,
			     'simplemin', 1))
      self.registerPMObject("NoseNVTLeapfrog",
                       ('temp', 300, 
                        'inertia', 0.5,
                        'bathpos', 1.0))
      self.registerPMObject("NPTVerlet",
                       ('temperature', 300,
		        'pressure', 1.0,
			'tauT', 1.0,
			'tauV', 1.0,
			'tauP', 1.0))
      self.registerPMObject("NumericalDifferentiation",
			    ('epsilon', 1))
      self.registerPMObject("PLeapfrog", ())
      self.registerPMObject("DMDLeapfrog",
                       ('iter', 5,
                        'gamma', 0.5,
                        'temp', 300,
                        'seed', 1234
                       ))
      self.registerPMObject("OpenMM",
                       ('temp', 300,
      		        'gamma', 2,
      			'seed', 1234,
      			'HarmonicBondForce', True,
      			'HarmonicAngleForce', True,
      			'RBDihedralForce', True,
      			'PeriodicTorsion', True,
      			'NonbondedForce', True,
      			'IntegratorType', 1,
      			'GBSAEpsilon', 1.0,
                        'GBSASolvent', 78.3,
                        'commonmotion', 0))

   


propFactory = PropagatorFactory()    #: PropagatorFactory singleton object

#objects = []
