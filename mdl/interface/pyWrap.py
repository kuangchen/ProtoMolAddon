import sys
import os
#from makeInterface import *


# List of Python modules which should be generated.
# The mapping is from a directory in the ProtoMol framework
# to a list of module names.
# The module names will correspond to prefixes for the interface
# file, source file and shared object.
mdlmodules = {'integrator/leapfrog':['LeapfrogIntegrator', 'LeapfrogTruncatedShadow', 'PLeapfrogIntegrator', 'DMDLeapfrogIntegrator', 'NoseNVTLeapfrogIntegrator'],
              'integrator/':['PySTSIntegrator', 'PyMTSIntegrator'],
	      'integrator/base':['LangevinImpulseIntegrator', 'CGMinimizerIntegrator', 'NumericalDifferentiation', 'LangevinLeapfrogIntegrator'],
              'type':['Vector3DBlock', 'ScalarStructure'],
	      'integrator/normal':['NormalModeBrownian', 'NormalModeDiagonalize', 'NormalModeMinimizer', 'NormalModeLangevin', 'NormalModeUtilities', 'NormalModeMori', 'NormalModeRelax'],
              'integrator/hessian':['HessianInt'],
              'io':['DCDTrajectoryReader', 'EigenvectorReader', 'EigenvectorTextReader', 'PARReader', 'PDBReader', 'PDBWriter', 'PSFReader', 'XYZBinReader', 'XYZReader', 'XYZTrajectoryReader', 'XYZTrajectoryWriter', 'XYZWriter'],
              'modifier':['ModifierShake'],
              'output':['OutputCache', 'OutputDCDTrajectory', 'OutputDCDTrajectoryVel', 'OutputEnergies', 'OutputFAHGUI', 'OutputFinalPDBPos', 'OutputFinalXYZPos', 'OutputFinalXYZVel', 'OutputScreen', 'OutputXYZTrajectoryForce', 'OutputXYZTrajectoryPos', 'OutputXYZTrajectoryVel'],
              'base':['MathUtilities'],
	      'topology':['TopologyUtilities', 'GenericTopology'],
	      'force/bonded':['BondForce','AngleForce','DihedralForce','HarmDihedralForce','ImproperForce'],
              'force/nonbonded':['SimpleFullForce', 'CutoffForce'],
              'force/system':['PySystemForce'],
              'force':['ForceGroup'],
              './force/nonbonded':['EwaldForce'],
              './':['ProtoMolApp']
	     }


def pyWrap(env, pmhome):
  import distutils.sysconfig
  env.Append(SWIGFLAGS=['-c++', '-python', '-w312', '-w314', '-w315', '-w317', '-w361', '-w362', '-w389', '-w401', '-w454', '-w503', '-w509', '-I'+pmhome],
                     CPPPATH=[distutils.sysconfig.get_python_inc()],
                     SHLIBPREFIX="",
                     ENV={'PATH':os.environ['PATH']})

  pyvers = "python%d.%d" % sys.version_info[:2]
  numpypath = sys.exec_prefix+"/lib/"+pyvers+"/site-packages/numpy/core/include/numpy/"

  env['ARFLAGS'] = "rcS"
  env.Append(CPPPATH = '#')
  env.Append(LINKFLAGS=' -Wl,-E')
  env.Append(CXXFLAGS=' -fPIC')
  env.Append(CXXFLAGs='-I'+pmhome)
  env.Append(SHCXXFLAGS=' -I'+numpypath)
  env.Append(SHCXXFLAGS='-I'+pmhome)
  env.Append(SHLINKFLAGS=' -Wl,-E')
  env.Append(_LIBDIRFLAGS="-L.")

  if (str(env['CCFLAGS']).find('HAVE_GUI') != -1):
    env.Append(SWIGFLAGS='-DHAVE_GUI')
  for dir in mdlmodules.iterkeys():
    modulelist = mdlmodules[dir]
    for i in range(0, len(modulelist)):
       module = modulelist[i]
       env.SharedLibrary(target=dir+'/_'+module+'.so', source=[dir+'/'+module+'.i'], SHLIBPREFIX="")
