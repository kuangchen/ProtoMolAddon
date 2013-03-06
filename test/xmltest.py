from xml.dom.minidom import parseString
import elementtree.ElementTree as ET
import glob
import os
import subprocess
import comparator
import argparse

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

# Arguments
parser = argparse.ArgumentParser(description='ProtoMol Test Suite')
parser.add_argument('--verbose', '-v', action='store_true', default=False, help='Verbose output')
parser.add_argument('--parallel', '-p', action='store_true', default=False, help='MPI Testing')

args = parser.parse_args()

# Setup Results XML
xml_root = ET.Element("TestRun")
xml_fail = ET.SubElement(xml_root, "FailedTests")
xml_pass = ET.SubElement(xml_root, "SuccessfulTests")

# Setup Statistics
stats_test = 0
stats_pass = 0
stats_fail = 0
stats_error = 0

# Find Executable
executable = []
if args.parallel:
    executable.append( "mpirun" )
    executable.append( "-np" )
    executable.append( "2" )
executable.append( os.path.join(os.getcwd(), 'ProtoMol') )

# Find Tests
tests = glob.glob(os.path.join(os.getcwd(), 'tests', '*.conf'))
tests.sort()

stats_test = len(tests)

# Run Tests
testid = 0
for test in tests:
    conf_param_overrides = parse_params(test)
    epsilon = conf_param_overrides.get('epsilon', DEFAULT_EPSILON)
    scaling_factor = conf_param_overrides.get('scaling_factor',
        DEFAULT_SCALINGFACTOR)

    testid += 1
    name = os.path.splitext(os.path.basename(test))[0]

    # Run ProtoMol
    cmd = executable[:]
    cmd.append( test )

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()

    if p.returncode > 0:
        stats_error += 1
        continue

    # Find Outputs
    expects = []
    outputs = glob.glob('tests/output/' + name + '.*')
    outputtemp = []
    for output in outputs:
        outbase = os.path.basename(output)
        if os.path.exists('tests/expected/' + outbase):
            outputtemp.append(output)
            expects.append('tests/expected/' + outbase)

    outputs = outputtemp[:]

    # Remove Untestable Types
    def untestable(path):
        ftype = os.path.splitext(os.path.basename(path))[1]
        if ftype in ['.dcd', '.header', '.xtc']:
            return False
        else:
            return True

    expects = filter(untestable, expects)
    outputs = filter(untestable, outputs)

    # Compare Outputs
    output_size = len(outputs)
    output_pass = []
    output_fail = []

    for i in xrange(output_size):
        ignoreSign = False
        ftype = os.path.splitext(os.path.basename(outputs[i]))[1]

        if ftype == '.vec':
            ignoreSign = True

        if comparator.compare(expects[i], outputs[i], epsilon, scaling_factor, ignoreSign):
            string = outputs[i] + " matches"
            output_pass.append( string )
            print string
        else:
            string = outputs[i] + " differs"
            output_fail.append( string )
            print string

    print name + " " + str(output_size) + " " + str(len(output_pass))

    # Create XML
    if len(output_pass) == output_size:
        stats_pass += 1
        xml_test = ET.SubElement(xml_pass,"Test")
        xml_test.set("id", str(testid))

        xml_test_name = ET.SubElement(xml_test,"Name")
        xml_test_name.text = name
    else:
        stats_fail += 1
        xml_test = ET.SubElement(xml_fail,"FailedTest")
        xml_test.set("id", str(testid))

        xml_test_name = ET.SubElement(xml_test,"Name")
        xml_test_name.text = name

        xml_test_type = ET.SubElement(xml_test,"FailureType")
        xml_test_type.text = "Assertion"

        def combine(a,b):
            return a + "\n" + b

        xml_test_message = ET.SubElement(xml_test,"Message")
        xml_test_message.text = reduce(combine, output_fail)

# Create Test Statistics
xml_stat = ET.SubElement(xml_root, "Statistics")

xml_stat_test = ET.SubElement(xml_stat, "Tests")
xml_stat_test.text = str(stats_test)

xml_stat_failt = ET.SubElement(xml_stat, "FailuresTotal")
xml_stat_failt.text = str(stats_fail)

xml_stat_error = ET.SubElement(xml_stat, "Errors")
xml_stat_error.text = str(stats_error)

xml_stat_fail = ET.SubElement(xml_stat, "Failures")
xml_stat_fail.text = str(stats_fail)

# Write Results
results = ET.tostring(xml_root)

fResult = open("results.xml", "w")
fResult.write( parseString(results).toprettyxml() )