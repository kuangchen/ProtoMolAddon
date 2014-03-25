#include <sstream>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/tof/OutputCEMRecorder.h>
#include <protomol/force/bonded/MTorsionSystemForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/parallel/Parallel.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/Report.h>

#include <string>
#include <algorithm>
#include <iterator>

using namespace ProtoMol;
using namespace ProtoMolAddon::Constant;
using namespace ProtoMolAddon::ToF;

using std::string;
using std::cerr;
using std::ifstream;


const string OutputCEMRecorder::keyword("OutputCEMRecorder");

OutputCEMRecorder::OutputCEMRecorder(): Output(-1) {}

OutputCEMRecorder::~OutputCEMRecorder() {}

OutputCEMRecorder::OutputCEMRecorder(const string &output_filename, double detector_pos) : 
  Output(1), recorder(detector_pos), output_filename(output_filename) { 
  cout << "detector_pos = " << detector_pos << "\n"; }

Output* OutputCEMRecorder::doMake(const std::vector<Value> &values) const
{
  return new OutputCEMRecorder(values[0], values[1]);
}

void OutputCEMRecorder::getParameters(std::vector<Parameter> &parameter) const
{
  parameter.push_back(Parameter(getId(), Value(output_filename, ConstraintValueType::NotEmpty())));
  //Output::getParameters(parameter);
  parameter.push_back(Parameter("DetectorPos", Value(detector_pos), ConstraintValueType::NoConstraints()));//;, 1.0,
//				Text("Detector position")));
}

void OutputCEMRecorder::doInitialize() {
  recorder.Initialize(app);
}

void OutputCEMRecorder::doRun(int step) {
  recorder.UpdateRecord(app);
}


void OutputCEMRecorder::doFinalize(int step) {
  try {
    ofstream f(output_filename);
    f << recorder;
  }
  catch (ofstream::failure e) {
    std::cerr << " Exception opening/reading/closing file " << output_filename << "\n";
  }
}
