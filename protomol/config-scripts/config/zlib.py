from SCons.Script import *
import config


def configure(conf):
    env = conf.env
    config.check_home(conf, 'zlib')
    config.check_home(conf, 'zlib', lib_suffix = '', inc_suffix = '/src')

    return config.check_header(conf, 'zlib.h') and config.check_lib(conf, 'z')
