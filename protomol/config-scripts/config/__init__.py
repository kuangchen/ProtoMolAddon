# Child modules
import compiler
import boost
import python
import xml
import pthreads
import expat
import openssl
import dist
import resources
import build_info
import packager
import nsi
import cx_Freeze
import run_distutils
import dmg
import valgrind
import malloc
import freetype2
import opengl
import rpm
import deb
import mkl
import osx
import lapack
import libsqlite
import pkg
import distpkg
import app
import gromacs
import glew
import glut
import zlib
import libbzip2

import sys
import os
import imp
import traceback
from SCons.Script import *


def check_env(name, require = False):
    if os.environ.has_key(name): return os.environ[name]
    if require: raise Exception, "Missing environment variable: " + name
    

def check_home(conf, name, inc_suffix = '/include', lib_suffix = '/lib',
               require = False, suffix = '_HOME'):
    name = name.upper()

    # Include
    inc = check_env(name + '_INCLUDE')
    if inc:
        conf.env.AppendUnique(CPPPATH = [inc])
        require = False

    # Lib path
    libpath = check_env(name + '_LIBPATH')
    if libpath:
        conf.env.Prepend(LIBPATH = [libpath])
        require = False

    # Link flags
    linkflags = check_env(name + '_LINKFLAGS')
    if linkflags:
        conf.env.AppendUnique(LINKFLAGS = linkflags.split())

    # Home
    home = check_env(name + suffix, require)
    if home:
        if not inc: conf.env.AppendUnique(CPPPATH = [home + inc_suffix])
        if not libpath: conf.env.Prepend(LIBPATH = [home + lib_suffix])


def check_lib(conf, lib, unique = False, append = False, **kwargs):
    if conf.CheckLib(lib, autoadd = 0, **kwargs):
        if unique:
            conf.env.Prepend(LIBS = [lib])
        else:
            conf.env.Prepend(LIBS = [lib])

        return True

    return False


def require_lib(conf, lib, **kwargs):
    if not check_lib(conf, lib, **kwargs):
        raise Exception, 'Need library ' + lib


def check_header(conf, hdr, **kwargs):
    return conf.CheckHeader(hdr, **kwargs)


def require_header(conf, hdr, **kwargs):
    if not check_header(conf, hdr, **kwargs):
        raise Exception, 'Need header ' + hdr


def check_cxx_header(conf, hdr, **kwargs):
    return conf.CheckCXXHeader(hdr, **kwargs)


def require_cxx_header(conf, hdr, **kwargs):
    if not check_cxx_header(conf, hdr, **kwargs):
        raise Exception, 'Need C++ header ' + hdr


def load_conf_module(name, home = None):
    modname = 'config.' + name;

    if not modname in sys.modules:
        env_name = name.upper().replace('-', '_')
        env_home = env_name + '_HOME'

        home = check_env(env_home)
        if home and not os.path.isdir(home):
            raise Exception, '$%s=%s is not a directory' % (env_home, home)

        filename = name + '-config.py'
        if home: filename = home + '/' + filename

        if os.path.exists(filename):
            try:
                imp.load_source(modname, filename)

                if home:
                    mod = sys.modules[modname]
                    mod.__dict__['home'] = home

            except Exception, e:
                raise Exception, "Error loading configuration file '%s' " \
                    "for '%s':\n%s" % (filename, name, e)
        else:
            raise Exception, "Configuration file '%s' not found " \
                "for '%s' not found please set %s" % (filename, name, env_home)

    return sys.modules[modname]



def call_single(mod, func, required = False, args = [], kwargs = {}):
    if type(mod) == str:
        try:
            mod = load_conf_module(mod)
        except:
            if required: raise
            return False

    try:
        function = mod.__dict__[func]
    except KeyError:
        if required:
            raise Exception, "function '%s' not found in '%s'" % (
                func, mod.__name__)

        return False

    try:
        return function(*args, **kwargs)

    except Exception, e:
        if required: raise
        return False


def call_conf_module(mods, func, required = False, args = [], kwargs = {}):
    ret = True

    if not isinstance(mods, list): mods = [mods]
    for mod in mods:
        if not call_single(mod, func, required, args, kwargs):
            ret = False

    return ret


def get_deps(mods):
    if not isinstance(mods, list): mods = [mods]

    deps = set()
    for mod in mods:
        if not mod in deps:
            deps.add(mod)

            try:
                mod = load_conf_module(mod)
                mods += mod.__dict__['deps']

            except KeyError:
                pass

            except Exception:
                pass

    return list(deps)


def add_vars(mods, vars, required = False):
    if not isinstance(mods, list): mods = [mods]
    mods = get_deps(mods)

    return call_conf_module(mods, 'add_vars', required, [vars])


def configure(mods, conf, required = True, **kwargs):
    return call_conf_module(mods, 'configure', required, [conf], kwargs)


def make_env(deps, extra_vars = None):
    configs = []

    if os.environ.has_key('SCONS_OPTIONS'):
        options = os.environ['SCONS_OPTIONS']
        if not os.path.exists(options):
            print 'options file "%s" set in SCONS_OPTIONS does not exist' % \
                options
            Exit(1)

        configs.append(options)

    if os.path.exists('default-options.py'):
        configs.append('default-options.py')
    if os.path.exists('options.py'): configs.append('options.py')

    vars = Variables(configs)
    add_vars(deps, vars)
    if extra_vars:
        for var in extra_vars: vars.AddVariables(var)
    env = Environment(variables = vars, ENV = os.environ)
    Help(vars.GenerateHelpText(env))

    Export('env')

    return env

