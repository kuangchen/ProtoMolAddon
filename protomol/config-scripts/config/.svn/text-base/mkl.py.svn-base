import sys
from SCons.Script import *
import config
from platform import machine, architecture


def CheckMKL(context):
    env = context.env
    conf = context.sconf

    context.Message('Checking for Intel MKL... ')

    # Save
    save_LIBPATH = env.get('LIBPATH', None)
    save_CPPPATH = env.get('CPPPATH', None)
    save_LIBS = env.get('LIBS', None)

    config.check_home(conf, 'mkl', suffix = 'ROOT')

    # Test program
    source = """
      #include "mkl.h"
      int main()
      {
        char ta = 'N', tb = 'N';
        int M = 1, N = 1, K = 1, lda = 1, ldb = 1, ldc = 1;
        double alpha = 1.0, beta = 1.0, *A = 0, *B = 0, *C = 0;
        dgemm(&ta, &tb, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
        return 0;
      }
    """
    source = source.strip()

    env.PrependUnique(LIBS = ['m'])

    # Compile & Link options
    compiler = env.get('compiler')
    if compiler == 'intel':
        if env['PLATFORM'] == 'win32':
            env.AppendUnique(CCFLAGS = ['/Qmkl'])
            env.AppendUnique(LINKFLAGS = ['/Qmkl'])
        else:
            env.AppendUnique(CCFLAGS = ['-mkl'])
            env.AppendUnique(LINKFLAGS = ['-mkl'])

    else:
        if architecture()[0] == '64bit': suffix = '_lp64'
        elif env['PLATFORM'] == 'win32': suffix = '_c'
        else: suffix = ''
        libs = ['mkl_intel' + suffix]

        if compiler == 'gnu': libs += ['mkl_gnu_thread', 'pthread']

        if env['PLATFORM'] == 'win32': libs += ['libiomp5mt']
        else: libs += ['iomp5']

        libs += ['mkl_core']

        env.Prepend(LIBS = libs)
        env.Prepend(LIBS = libs) # Add twice

    # Try it
    if context.TryCompile(source, '.cpp') and context.TryLink(source, '.cpp'):
        env.AppendUnique(CPPDEFINES = ['HAVE_MKL'])
        context.Result(True)
        return True

    # Restore
    env.Replace(LIBPATH = save_LIBPATH)
    env.Replace(CPPPATH = save_CPPPATH)
    env.Replace(LIBS = save_LIBS)

    context.Result(False)
    return False


def configure(conf):
    conf.AddTest('CheckMKL', CheckMKL)
    return conf.CheckMKL()
