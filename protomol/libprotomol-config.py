import os
from SCons.Script import *
import config
from platform import machine, architecture

deps = ['libfah', 'mkl', 'lapack', 'gromacs']


def add_vars(vars):
    vars.AddVariables(
        BoolVariable('fah', 'Set to 1 to build library for Folding@home', 0),
        BoolVariable('qrdiag', 'Set to 1 if QR diagonalization', 0),
        BoolVariable('gui', 'Set to 1 if using the GUI', 0),
        EnumVariable('lapack', 'Use LAPACK', 'any', allowed_values =
                     ('1', 'any', 'none', 'mkl', 'simtk', 'system')),
        BoolVariable('openmm', 'Build with OpenMM support', 0),
        BoolVariable('ltmdopenmm', 'Build with LTMD OpenMM support', 0),
        BoolVariable('gromacs', 'Enable or disable gromacs support', 0),
        BoolVariable('gromacs_standard', 'Enable or disable gromacs support',
                     0),
        )


def configure_deps(conf):
    env = conf.env

    # libfah
    if env.get('fah', 0): config.configure('libfah', conf)

    # DIAG Options
    if env.get('qrdiag', 0): env.AppendUnique(CPPDEFINES = ['HAVE_QRDIAG'])

    # GUI Options
    if env.get('gui',0):
        if env['PLATFORM'] == 'win32': config.require_lib(conf, 'wsock32')
        else: config.require_lib(conf, 'pthread')

        env.AppendUnique(CPPDEFINES = ['HAVE_GUI'])

    # LAPACK
    have_lapack = False
    lapack = env.get('lapack', 'any')
    if lapack == '1' or lapack is True: lapack = 'any'
    elif lapack == '0' or lapack is False: lapack = 'none'

    if lapack != 'none':
        # Intel MKL LAPACK
        if not have_lapack and lapack in ['any', 'mkl']:
            have_lapack = config.configure('mkl', conf)
            if have_lapack:
                env.AppendUnique(CPPDEFINES = ['HAVE_MKL_LAPACK'])
            elif lapack == 'mkl':
                raise Exception, "Missing MKL LAPACK"

        # System LAPACK
        if not have_lapack and lapack in ['any', 'system']:
            have_lapack = config.configure('lapack', conf)
            if not have_lapack and lapack == 'lapack':
                raise Exception, "Missing LAPACK"

        # SimTK LAPACK
        if not have_lapack and lapack in ['any', 'simtk']:
            config.check_home(conf, 'simtk_lapack')

            if (config.check_lib(conf, 'SimTKlapack') and
                config.check_cxx_header(conf, 'SimTKlapack.h')):

                env.AppendUnique(CPPDEFINES = ['HAVE_SIMTK_LAPACK'])
                have_lapack = True
 
            elif lapack == 'simtk_lapack':
                raise Exception, "Missing SimTK LAPACK"

        if not have_lapack: raise Exception, "Missing LAPACK"
						
    # OpenMM
    openmm = env.get('openmm', 0)
    if openmm:
        home = config.check_env('OPENMM_HOME', True)

        conf.env.AppendUnique(CPPPATH = [home + 'olla/include'])
        conf.env.AppendUnique(CPPPATH = [home + 'openmmapi/include'])
        config.require_cxx_header(conf, 'OpenMM.h')
        config.require_cxx_header(conf, 'openmm/Kernel.h')

        conf.env.Prepend(LIBPATH = [home])
        config.require_lib(conf, 'OpenMM')

        env.AppendUnique(CPPDEFINES = ['HAVE_OPENMM'])
						
    # LTMD OpenMM
    ltmd = env.get('ltmdopenmm', 0)
    if ltmd and openmm:
        config.check_home(conf, 'ltmdopenmm')
        config.require_lib(conf, 'OpenMMLTMD')
        env.AppendUnique(CPPDEFINES = ['HAVE_OPENMM_LTMD'])
			

    # Gromacs
    gromacs = env.get('gromacs', 0)
    if gromacs: config.configure('gromacs', conf)


    # Gromacs Standard
    gromacs_standard = env.get('gromacs_standard', 0)
    if gromacs_standard:
        config.check_home(conf, 'gromacs')
        config.require_lib(conf, 'md')
        config.require_lib(conf, 'gmx')
        env.AppendUnique(CPPDEFINES = ['HAVE_GROMACS'])


def configure(conf):
    env = conf.env

    configure_deps(conf)

    # Libprotomol paths
    if home:
        env.Append(CPPPATH = [home + '/src'])
        env.Prepend(LIBPATH = [home])

    # Library name
    lib = 'protomol'
    if env.get('fah', 0): lib += '-fah'

    config.require_lib(conf, lib)
