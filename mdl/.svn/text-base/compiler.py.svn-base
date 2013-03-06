import copy
import re
import os
from platform import machine
from SCons.Script import *
from subprocess import *
from SCons.Util import MD5signature

class static_lib_decider_hack:
    def __init__(self, env):
        self.env = env
        self.decider = env.decide_source

    def __call__(self, dep, target, prev_ni):
        if str(dep).endswith(self.env['LIBSUFFIX']):
            try:
                csig = dep.csig
            except AttributeError:
                csig = MD5signature(dep.get_contents())
                dep.csig = csig
                dep.get_ninfo().csig = csig

            #print dependency, csig, "?=", prev_ni.csig
            if prev_ni is None: return True
            return csig != prev_ni.csig

        return self.decider(dep, target, prev_ni)


def add_vars(vars):
    vars.AddVariables(
        ('optimize', 'Enable or disable optimizations', -1),
        ('globalopt', 'Enable or disable global optimizations', 0),
        ('sse2', 'Enable SSE2 instructions', 0),
        ('sse3', 'Enable SSE3 instructions', 0),
        ('auto_dispatch', 'Enable auto-dispatch of optimized code paths', 1),
        BoolVariable('debug', 'Enable or disable debug options',
                     os.getenv('DEBUG_MODE', 0)),
        BoolVariable('strict', 'Enable or disable strict options', 1),
        BoolVariable('threaded', 'Enable or disable thread support', 1),
        BoolVariable('profile', 'Enable or disable profiler', 0),
        BoolVariable('depends', 'Enable or disable dependency files', 0),
        BoolVariable('distcc', 'Enable or disable distributed builds', 0),
        BoolVariable('ccache', 'Enable or disable cached builds', 0),
        EnumVariable('platform', 'Override default platform', '',
                   allowed_values = ('', 'win32', 'posix')),
        EnumVariable('cxxstd', 'Set C++ language standard', 'gnu++98',
                   allowed_values = ('gnu++98', 'c++98', 'c++0x')),
        EnumVariable('compiler', 'Select compiler', 'default',
                   allowed_values = ('default', 'gnu', 'intel', 'mingw', 'msvc',
                                     'linux-mingw', 'aix', 'posix', 'hp', 'sgi',
                                     'sun')),
        BoolVariable('static', 'Link to static libraries', 0),
        BoolVariable('mostly_static', 'Link most libraries statically', 0),
        ('num_jobs', 'Set the concurrency level.', -1),
        )


def configure(conf, c99_mode = 1):
    env = conf.env

    if env.GetOption('clean'): return

    # Static lib decider hack
    env.Decider(static_lib_decider_hack(env))

    # Get options
    debug = env.get('debug')
    optimize = env.get('optimize')
    if optimize == -1: optimize = not debug
    globalopt = env.get('globalopt')
    sse2 = int(env.get('sse2', 0))
    sse3 = int(env.get('sse3', 0))
    auto_dispatch = int(env.get('auto_dispatch', 1))
    strict = int(env.get('strict', 1))
    threaded = int(env.get('threaded', 1))
    profile = int(env.get('profile', 0))
    depends = int(env.get('depends', 0))
    compiler = env.get('compiler', 0)
    distcc = env.get('distcc', 0)
    ccache = env.get('ccache', 0)
    cxxstd = env.get('cxxstd', 'c++0x')
    platform = env.get('platform', '')
    static = int(env.get('static', 0))
    mostly_static = int(env.get('mostly_static', 0))
    num_jobs = env.get('num_jobs', -1)

    if platform != '': env.Replace(PLATFORM = platform)

    # Select compiler
    compiler_mode = None

    # Prefer Intel compiler
    if compiler == 'default' and os.environ.get('INTEL_LICENSE_FILE', False):
        compiler = 'intel'

    if compiler:
        if compiler == 'gnu':
            Tool('gcc')(env)
            Tool('g++')(env)
            compiler_mode = "gnu"

        elif compiler == 'intel':
            Tool('intelc')(env)
            env['ENV']['INTEL_LICENSE_FILE'] = (
                os.environ.get("INTEL_LICENSE_FILE", ''))

            if env['PLATFORM'] == 'win32': compiler_mode = 'msvc'
            else: compiler_mode = 'gnu'

            if compiler_mode == 'msvc':
                env.Replace(AR = 'xilib')

                # Work around double CCFLAGS bug
                env.Replace(CXXFLAGS = ['/TP'])

        elif compiler == 'linux-mingw':
            env.Replace(CC = 'i586-mingw32msvc-gcc')
            env.Replace(CXX = 'i586-mingw32msvc-g++')
            env.Replace(RANLIB = 'i586-mingw32msvc-ranlib')
            env.Replace(PROGSUFFIX = '.exe')
            compiler_mode = "gnu"

        elif compiler == 'posix':
            Tool('cc')(env)
            Tool('cxx')(env)
            Tool('link')(env)
            Tool('ar')(env)
            Tool('as')(env)
            compiler_mode = "unknown"

        elif compiler in ['hp', 'sgi', 'sun', 'aix']:
            Tool(compiler + 'cc')(env)
            Tool(compiler + 'c++')(env)
            Tool(compiler + 'link')(env)

            if compiler in ['sgi', 'sun']:
                Tool(compiler + 'ar')(env)

            compiler_mode = "unknown"

        elif compiler != 'default':
            Tool(compiler)(env)


    if compiler_mode is None:
        if env['CC'] == 'cl' or env['CC'] == 'icl': compiler_mode = 'msvc'
        elif env['CC'] == 'gcc' or env['CC'] == 'icc': compiler_mode = 'gnu'
        else: compiler_mode = 'unknown'

    if compiler == 'default':
        cc = env['CC']
        if cc == 'cl': compiler = 'msvc'
        elif cc == 'gcc': compiler = 'gnu'
        elif cc == 'icl' or cc == 'icc': compiler = 'intel'


    print "Compiler: " + env['CC'] + ' (%s)' % compiler
    print "Platform: " + env['PLATFORM']
    print "Mode: " + compiler_mode


    # Options
    if compiler_mode == 'msvc':
        env.Append(CCFLAGS = ['/EHsc'])
        env.Append(CCFLAGS = ['/wd4297', '/wd4103'])
        env.Append(CPPDEFINES = ['_CRT_SECURE_NO_WARNINGS'])

        if compiler == 'intel':
            env.Append(CCFLAGS = ['/wd1786'])


    # Profiler flags
    if profile:
        if compiler_mode == 'gnu':
            env.Append(CCFLAGS = ['-pg'])
            env.Append(LINKFLAGS = ['-pg'])


    # Debug flags
    if debug:
        if compiler_mode == 'msvc':
            env.Append(CCFLAGS = ['/W1', '/Zi'])
            env.Append(LINKFLAGS = ['/DEBUG', '/MAP:${TARGET}.map'])
            env['PDB'] = '${TARGET}.pdb'

        elif compiler_mode == 'gnu':
            if compiler == 'gnu':
                env.Append(CCFLAGS = ['-ggdb', '-Wall'])
                env.Append(LINKFLAGS = ['-rdynamic']) # for backtrace
            elif compiler == 'intel':
                env.Append(CCFLAGS = ['-g', '-diag-enable', 'warn'])

            if strict: env.Append(CCFLAGS = ['-Werror'])

        env.Append(CPPDEFINES = ['DEBUG'])

    else:
        if compiler_mode == 'gnu':
            # Strip symbols
            env.Append(LINKFLAGS = ['-Wl,-s'])

        env.Append(CPPDEFINES = ['NDEBUG'])


    # Optimizations
    if optimize:
        # Machine
        if machine() != 'x86_64' and not (sse2 or sse3):
            if compiler == 'intel':
                if compiler_mode == 'gnu':
                    env.Append(CCFLAGS = ['-mia32'])
                elif compiler_mode == 'msvc':
                    env.Append(CCFLAGS = ['/arch:IA32'])
            elif compiler == 'gnu':
                env.Append(CCFLAGS = ['-march=pentium3'])
            elif compiler == 'msvc':
                env.Append(CCFLAGS = ['/arch:SSE'])

        if compiler_mode == 'gnu':
            env.Append(CCFLAGS = ['-O3', '-funroll-loops'])
            if compiler != 'intel':
                env.Append(CCFLAGS = ['-mfpmath=sse', '-ffast-math',
                                      '-fno-unsafe-math-optimizations'])
        elif compiler_mode == 'msvc':
            env.Append(CCFLAGS = ['/Ox'])
            if compiler == 'intel' and not globalopt:
                env.Append(LINKFLAGS = ['-qnoipo'])

        # Whole program optimizations
        if globalopt:
            if compiler == 'intel':
                if compiler_mode == 'gnu':
                    env.Append(LINKFLAGS = ['-ipo'])
                    env.Append(CCFLAGS = ['-ipo'])
                elif compiler_mode == 'msvc':
                    env.Append(LINKFLAGS = ['/Qipo'])
                    env.Append(CCFLAGS = ['/Qipo'])

            elif compiler == 'msvc':
                env.Append(CCFLAGS = ['/GL'])
                env.Append(LINKFLAGS = ['/LTCG'])
                env.Append(ARFLAGS = ['/LTCG'])

        # SSE optimizations
        if sse2:
            if compiler_mode == 'gnu':
                if compiler == 'intel':
                    env.Append(CCFLAGS = ['-xsse2']);
                    if auto_dispatch:
                        env.Append(CCFLAGS = ['-axSSE3,SSSE3,SSE4.1,SSE4.2'])
                else:
                    env.Append(CCFLAGS = ['-msse2']);
            elif compiler_mode == 'msvc':
                if compiler == 'intel':
                    env.Append(CCFLAGS = ['/QxSSE2']);
                    if auto_dispatch:
                        env.Append(CCFLAGS = ['/QaxSSE3,SSSE3,SSE4.1,SSE4.2'])
                else:
                    env.Append(CCFLAGS = ['/arch:SSE2']);

        elif sse3:
            if compiler_mode == 'gnu':
                if compiler == 'intel':
                    env.Append(CCFLAGS = ['-xsse3']);
                    if auto_dispatch:
                        env.Append(CCFLAGS = ['-axSSSE3,SSE4.1,SSE4.2'])
                else:
                    env.Append(CCFLAGS = ['-msse3']);
            elif compiler_mode == 'msvc':
                if compiler == 'intel':
                    env.Append(CCFLAGS = ['/QxSSE3']);
                    if auto_dispatch:
                        env.Append(CCFLAGS = ['/QaxSSSE3,SSE4.1,SSE4.2'])
                else:
                    env.Append(CCFLAGS = ['/arch:SSE3']);

        elif compiler == 'intel' and auto_dispatch:
            if compiler_mode == 'gnu':
                env.Append(CCFLAGS = ['-axSSE2,SSE3,SSSE3,SSE4.1,SSE4.2'])
            elif compiler_mode == 'msvc':
                env.Append(CCFLAGS = ['/QaxSSE2,SSE3,SSSE3,SSE4.1,SSE4.2'])


    # Pointer disambiguation
    if compiler == 'intel':
        if compiler_mode == 'gnu':
            env.Append(CCFLAGS = ['-restrict'])
        elif compiler_mode == 'msvc':
            env.Append(CCFLAGS = ['/Qrestrict'])


    # Dependency files
    if depends and compiler_mode == 'gnu':
        env.Append(CCFLAGS = ['-MMD -MF ${TARGET}.d'])


    # C mode
    if c99_mode:
        if compiler_mode == 'gnu':
            env.Append(CFLAGS = ['-std=c99'])
            env.Append(CXXFLAGS = ['-std=' + cxxstd])
        elif compiler_mode == 'msvc':
            env.Append(CFLAGS = ['/TP']) # C++ mode


    # Threads
    if threaded:
        if compiler_mode == 'gnu':
            if not conf.CheckLib('pthread'):
                print 'Need pthreads'
                Exit(1)

            env.Append(LINKFLAGS = ['-pthread'])
            env.Append(CPPDEFINES = ['_REENTRANT'])

        elif compiler_mode == 'msvc':
            if debug: env.Append(CCFLAGS = ['/MTd'])
            else: env.Append(CCFLAGS = ['/MT'])


    # Link flags
    if compiler_mode == 'msvc' and not optimize:
        env.Append(LINKFLAGS = ['/INCREMENTAL'])

    # static
    if static:
        if compiler_mode == 'gnu':
            env.Append(LINKFLAGS = ['-static'])


    # Num jobs default
    default_num_jobs = 1

    # distcc
    if distcc and compiler == 'gnu':
        default_num_jobs = 2
        env.Replace(CC = 'distcc ' + env['CC'])
        env.Replace(CXX = 'distcc ' + env['CXX'])

    # cccache
    if ccache and compiler == 'gnu':
        env.Replace(CC = 'ccache ' + env['CC'])
        env.Replace(CXX = 'ccache ' + env['CXX'])

    # Num jobs
    if num_jobs == -1:
        if os.environ.has_key('SCONS_JOBS'):
            num_jobs = int(os.environ.get('SCONS_JOBS', num_jobs))
        else: num_jobs = default_num_jobs

    SetOption('num_jobs', num_jobs)
    print "running with -j", GetOption('num_jobs')


    # For darwin
    if env['PLATFORM'] == 'darwin':
        env.Append(CPPDEFINES = ['__APPLE__'])
        env['PLATFORM'] = 'posix'


def findLibPath(env, lib):
    eenv = copy.copy(os.environ)
    eenv['LIBRARY_PATH'] = ':'.join(env['LIBPATH'])
    cmd = env['CXX'].split()
    libpat = env['LIBPREFIX'] + '%s' + env['LIBSUFFIX']

    path = Popen(cmd + ['-print-file-name=' + libpat % lib],
                 stdout = PIPE, env = eenv).communicate()[0].strip()

    if path == libpat % lib: return None
    return path


def mostly_static_libs(env, ignore = ['pthread', 'dl']):
    eenv = copy.copy(os.environ)
    eenv['LIBRARY_PATH'] = ':'.join(env['LIBPATH'])
    cmd = env['CXX'].split()
    libpat = env['LIBPREFIX'] + '%s' + env['LIBSUFFIX']

    libs = []

    for lib in env['LIBS']:
        skip = False
        for i in ignore:
            if re.match(i, lib):
                skip = True
                break

        if not (skip or lib.startswith(os.sep) or
                lib.endswith(env['LIBSUFFIX'])):
            path = Popen(cmd + ['-print-file-name=' + libpat % lib],
                         stdout = PIPE, env = eenv).communicate()[0].strip()
            if path == libpat % lib: libs.append(lib)
            else: libs.append(File(path))

        else: libs.append(lib)

    env.Replace(LIBS = libs)
