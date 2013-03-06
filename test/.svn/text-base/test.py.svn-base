#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import glob
import shlex
import subprocess
import comparator
import argparse

import logging

DEFAULT_EPSILON = 0.00001
DEFAULT_SCALINGFACTOR = 1.0


def parse_params(flname):
    """
....Parses test parameters from Protomol configuration files.
....Paramaters are used to change how testing is performed.
....It is assumed that parameters will be in the form:
            ## key = value
....where key will be a string, value will by any python structure
....that fits on one line, and the hashmarks will be the first two
        characters on the line.
...."""

    fl = open(flname)
    params = dict()
    for ln in fl:

        if ln[:2] != '##':
            continue
        eq_pos = ln.find('=')
        if eq_pos == -1:
            continue
        param_name = ln[2:eq_pos].strip()
        param_str = ln[eq_pos + 1:].strip()
        param_val = eval(param_str)
        params[param_name] = param_val
    fl.close()
    return params


def run_test(protomol_path, conf_file, pwd):
    tests = 0
    testspassed = 0
    testsfailed = 0
    failedtests = []

    conf_param_overrides = parse_params(conf_file)
    epsilon = conf_param_overrides.get('epsilon', DEFAULT_EPSILON)
    scaling_factor = conf_param_overrides.get('scaling_factor',
            DEFAULT_SCALINGFACTOR)

    base = os.path.splitext(os.path.basename(conf_file))[0]
    logging.info('Executing Test: ' + base)

    cmd = protomol_path[:]
    cmd.append(conf_file)

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    if p.returncode > 0:
        s = 'Not able to execute Protomol!\n'
        s += 'cmd: ' + str(cmd) + '\n'
        s += 'stdout: ' + stdout + '\n'
        s += 'stderr: ' + stderr + '\n'
        logging.critical(s)
    expects = []
    outputs = glob.glob('tests/output/' + base + '.*')
    outputtemp = []
    for output in outputs:
        outbase = os.path.basename(output)
        if os.path.exists('tests/expected/' + outbase):
            outputtemp.append(output)
            expects.append('tests/expected/' + outbase)

    outputs = outputtemp

    for i in xrange(len(outputs)):
        tests += 1
        ftype = os.path.splitext(os.path.basename(outputs[i]))[1]

        if ftype in ['.dcd', '.header', '.xtc']:
            continue

        ignoreSign = False

        # Ignore signs on eignevectors

        if ftype == '.vec':
            ignoreSign = True
        logging.info('Testing: ' + expects[i] + ' ' + outputs[i])

        if comparator.compare(expects[i], outputs[i], epsilon,
                              scaling_factor, ignoreSign):
            logging.info('Passed')
            testspassed += 1
        else:
            logging.warning('Failed')
            testsfailed += 1
            failedtests.append('Comparison of ' + expects[i] + ' and '
                               + outputs[i])

    return (tests, testspassed, testsfailed, failedtests)


def find_conf_files(pwd, args):
    files = []
    if args.single == None and args.regex == None:
        files = glob.glob(os.path.join(pwd, 'tests', '*.conf'))
    elif args.single != None:
        if args.single.find('.conf') != -1:
            files.append(args.single)
        else:
            raise Exception, args.single \
                + ' is not a valid configuration file'
    elif args.regex != None:
        files = glob.glob(os.path.join(pwd, 'tests', args.regex
                          + '.conf'))
    return files


def find_protomol(pwd):
    path = []
    if args.parallel:
        path.append('mpirun')
        path.append('-np')
        path.append('2')
    unix_path = os.path.join(pwd, 'ProtoMol')
    win_path = os.path.join(pwd, 'ProtoMol.exe')
    if os.path.exists(unix_path):
        path.append(unix_path)
        return path
    elif os.path.exists(win_path):
        path.append(win_path)
        return path
    else:
        raise Exception, 'Cannot find ProtoMol executable in ' + pwd \
            + ' . Please put the ProtoMol executable in this directory'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ProtoMol Test Suite')
    parser.add_argument('--verbose', '-v', action='store_true',
                        default=False, help='Verbose output')
    parser.add_argument('--parallel', '-p', action='store_true',
                        default=False, help='MPI Testing')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--single', '-s',
                       help='Single test to run. Must be within the tests directory.'
                       )
    group.add_argument('--regex', '-r',
                       help='Regular expression of tests to run, Requires quotation marks around argument'
                       )

    args = parser.parse_args()

    pwd = os.getcwd()

    level = logging.INFO
    if args.verbose:
        level = logging.DEBUG
    FORMAT = '%(message)s'
    logging.basicConfig(level=level, format=FORMAT)

    files = find_conf_files(pwd, args)
    files.sort()
    logging.debug('Files: ' + str(files))

    tests = 0
    testspassed = 0
    testsfailed = 0
    failedtests = []

    protomol_path = find_protomol(pwd)

    for conf_file in files:
        results = run_test(protomol_path, conf_file, pwd)
        tests += results[0]
        testspassed += results[1]
        testsfailed += results[2]
        failedtests.extend(results[3])

    testsnotrun = tests - (testspassed + testsfailed)

    logging.info('')
    logging.info('Tests: %d' % tests)
    logging.info('Tests Not Run: %d' % testsnotrun)
    logging.info('')
    logging.info('Tests Passed: %d' % testspassed)
    logging.info('Tests Failed: %d\n' % testsfailed)

    if len(failedtests) > 0:
        for test in failedtests:
            logging.warning('%s failed.' % test)
        sys.exit(1)
