if not 'PROTOMOL_HOME' in locals():
    PROTOMOL_HOME = '.'

execfile(PROTOMOL_HOME + '/../mdl/configfuncs.py')

def add_config_vars(vars):
    vars.AddVariables(
        BoolVariable('qrdiag', 'Set to 1 if QR diagonalization', 0),
        BoolVariable('gui', 'Set to 1 if using the GUI', 0),
        BoolVariable('lapack', 'Use LAPACK', 0),
        BoolVariable('simtk_lapack', 'Use SimTK LAPACK', 0),
        EnumVariable('openmm', 'Build with OpenMM', 'none',
                   allowed_values = ('none', 'reference', 'cuda'))
    )

def config_configure():
    # DIAG Options
    use_qr = int( env.get( 'qrdiag', 0 ) )
    if use_qr == 1:
        env.Append(CCFLAGS = '-DHAVE_QRDIAG')

    # GUI Options
    use_gui = int( env.get( 'gui',0 ) )
    if use_gui == 1:
        env.Append(CCFLAGS = '-DHAVE_GUI')

        if env['PLATFORM'] == 'win32':
            check_library( 'wsock32', True )
        else:
            check_library( 'pthread', True )

    # LAPACK
    have_lapack = False
    use_lapack = int( env.get( 'lapack', 0 ) )
    if use_lapack:
        lapack_home = check_envvar( 'LAPACK_HOME' )

        if lapack_home != None:
            env.Append(CPPPATH = [lapack_home])
            env.Append(LIBPATH = [lapack_home])

        if check_library( 'lapack' ):
            env.Append(CPPDEFINES = ['HAVE_LAPACK'])
            have_lapack = True

            if env['PLATFORM'] == 'posix':
                # BLAS
                home = check_envvar( 'BLAS_HOME' )
                if home != None: env.Append(LIBPATH = [home])
                check_library( 'blas' )

                # G2C
                home = check_envvar( 'G2C_HOME' )
                if home != None: env.Append(LIBPATH = [home])
                check_library( 'g2c' )

                # GFortran
                home = check_envvar( 'GFORTRAN_HOME' )
                if home != None: env.Append(LIBPATH = [home])
                check_library( 'gfortran' )

    # SimTK LAPACK
    use_simtk = int( env.get( 'simtk_lapack', 0) )
    if use_simtk == 1 or (use_lapack and not have_lapack):
        simtk_home = check_envvar( 'SIMTK_LAPACK_HOME' )

        if simtk_home != None:
            env.Append(LIBPATH = [simtk_home + '/lib'])
            env.Append(CPPPATH = [simtk_home + '/include'])

        if check_library( 'SimTKlapack' ) and \
                check_header( 'SimTKlapack.h', True ):
            env.Append(CPPDEFINES = ['HAVE_SIMTK_LAPACK'])
            have_lapack = True

    if not have_lapack and (use_simtk or use_lapack):
        print "Missing lapack"
        sys.exit(1)

    # OpenMM Options
    openmm_type = env.get( 'openmm', 'none' )
    if openmm_type != 'none':       
        # The following must bail if it is not found as openmm is not
        # installed to a place that the compiler will locate by default. The
        # same is also true for CUDA.
        openmm_home = check_envvar( 'OPENMM_HOME', True )

        env.Append(CPPPATH = [openmm_home + os.sep + 'include'])
        env.Append(LIBPATH = [openmm_home + os.sep + 'lib'    ])

        if openmm_type == 'reference':
            if check_library( 'OpenMM_d', True ):
                env.Append(CPPDEFINES = ['HAVE_OPENMM'])

        if openmm_type == 'cuda':
            cuda_home = check_envvar( 'CUDA_HOME', True )

            env.Append(LIBPATH = [cuda_home + os.sep + 'lib'])

            if check_library( 'OpenMM_d', True ) and \
                    check_library( 'OpenMMCuda_d', True ) and \
                    check_library( 'cudart', True ):
                env.Append(CPPDEFINES = ['HAVE_OPENMM'])
