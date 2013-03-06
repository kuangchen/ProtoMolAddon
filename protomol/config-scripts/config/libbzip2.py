from SCons.Script import *
import config


def configure(conf):
    env = conf.env
    config.check_home(conf, 'libbzip2')
    config.check_home(conf, 'libbzip2', lib_suffix = '', inc_suffix = '/src')

    return \
        config.check_header(conf, 'bzlib.h') and config.check_lib(conf, 'bz2')
