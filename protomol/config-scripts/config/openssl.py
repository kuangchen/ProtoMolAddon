import os
from SCons.Script import *
import config


def configure(conf):
    env = conf.env

    if os.environ.has_key('OPENSSL_HOME'):
        home = os.environ['OPENSSL_HOME']
        if os.path.exists(home + '/inc32') and os.path.exists(home + '/out32'):
            env.AppendUnique(CPPPATH = [home + '/inc32'])
            env.AppendUnique(LIBPATH = [home + '/out32'])
        else:
            env.AppendUnique(CPPPATH = [home + '/include'])
            if os.path.exists(home + '/lib'):
                env.AppendUnique(LIBPATH = [home + '/lib'])
            else:
                env.AppendUnique(LIBPATH = [home])

    if env['PLATFORM'] == 'posix': config.check_lib(conf, 'dl')

    if (conf.CheckCHeader('openssl/ssl.h') and
        (config.check_lib(conf, 'crypto') and
         config.check_lib(conf, 'ssl')) or
        (config.check_lib(conf, 'libeay32') and
         config.check_lib(conf, 'ssleay32'))):

        if env['PLATFORM'] == 'win32':
            for lib in ['wsock32', 'advapi32', 'gdi32', 'user32']:
                config.require_lib(conf, lib)
 
        env.AppendUnique(CPPDEFINES = ['HAVE_OPENSSL'])
        return True

    else: raise Exception, 'Need openssl'

