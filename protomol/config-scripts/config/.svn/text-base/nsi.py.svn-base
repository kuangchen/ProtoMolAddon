from SCons.Script import *
from SCons.Action import CommandAction


def build_function(target, source, env):
    nsi = str(source[0])

    env.Replace(package = str(target[0]))

    tmp = nsi + '.tmp'
    input = None
    output = None
    try:
        input = open(nsi, 'r')
        output = open(tmp, 'w')
        output.write(input.read() % env)

    finally:
        if input is not None: input.close()
        if output is not None: output.close()

    action = CommandAction('$NSISCOM $NSISOPTS ' + tmp)
    return action.execute(target, [tmp], env)


def configure(conf):
    env = conf.env

    env['NSISCOM'] = 'makensis.exe'
    env['NSISOPTS'] = ''

    bld = Builder(action = build_function)
    env.Append(BUILDERS = {'Nsis' : bld}) 
