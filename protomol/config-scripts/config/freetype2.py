import os
from SCons.Script import *
import config

deps = ['osx']

def configure(conf):
    env = conf.env

    if os.environ.has_key('FREETYPE2_INCLUDE'):
        env.AppendUnique(CPPPATH =
                         os.environ['FREETYPE2_INCLUDE'].split(os.pathsep))
    else:
        try:
            env.ParseConfig('freetype-config --cflags')
        except OSError: pass

    if os.environ.has_key('FREETYPE2_LIBPATH'):
        env.AppendUnique(LIBPATH = [os.environ['FREETYPE2_LIBPATH']])

    if env['PLATFORM'] == 'darwin':
        config.configure('osx', conf)
        if not conf.CheckOSXFramework('CoreServices'):
            raise Exception, 'Need CoreServices framework'
        if not conf.CheckOSXFramework('ApplicationServices'):
            raise Exception, 'Need ApplicationServices framework'

    config.require_header(conf, 'ft2build.h')
    config.require_lib(conf, 'freetype')

    config.check_home(conf, 'zlib', lib_suffix = '') # TODO Hack!
    config.require_lib(conf, 'z')

    return True
