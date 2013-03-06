import os
from SCons.Script import *
import config


def configure(conf):
    env = conf.env

    libname = 'expat'

    if os.environ.has_key('EXPAT_INCLUDE'):
        env.AppendUnique(CPPPATH = [os.environ['EXPAT_INCLUDE']])

    if os.environ.has_key('EXPAT_LIBPATH'):
        env.AppendUnique(LIBPATH = [os.environ['EXPAT_LIBPATH']])

    if env['PLATFORM'] == 'win32':
        libname = 'lib' + libname + 'MT'
        env.AppendUnique(CPPDEFINES = ['XML_STATIC'])

    if conf.CheckCHeader('expat.h') and config.check_lib(conf, libname):
        env.AppendUnique(CPPDEFINES = ['HAVE_EXPAT'])

    else: raise Exception, 'Need expat'

    return True
