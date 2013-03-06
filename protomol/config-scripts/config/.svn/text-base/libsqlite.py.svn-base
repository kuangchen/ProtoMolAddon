import os
import platform
from SCons.Script import *
import config

deps = ['pthreads']

def configure(conf):
    env = conf.env

    # pthread
    if env['PLATFORM'] != 'win32':
        config.configure('pthreads', conf)

    # wsock32
    if env['PLATFORM'] == 'win32':
        config.require_lib(conf, 'wsock32')

    config.check_home(conf, 'libsqlite', '', '')

    if config.check_lib(conf, 'sqlite3') and \
            config.check_header(conf, 'sqlite3.h'):
        env.AppendUnique(CPPDEFINES = ['HAVE_LIBSQLITE'])
        return True

    return False
