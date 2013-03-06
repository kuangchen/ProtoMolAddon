import os
from SCons.Script import *
import config

def configure(conf):
    env = conf.env

    if os.environ.has_key('PTHREADS_HOME'):
        home = os.environ['PTHREADS_HOME']
        env.AppendUnique(CPPPATH = [home])
        env.AppendUnique(LIBPATH = [home])

    if conf.CheckCHeader('pthread.h') and config.check_lib(conf, 'pthread'):
        env.AppendUnique(CPPDEFINES = ['HAVE_PTHREADS'])
    else: raise Exception, 'Need pthreads'
