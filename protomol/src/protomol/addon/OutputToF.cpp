#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <sstream>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/OutputToF.h>
#include <string>
#include <algorithm>
#include <iterator>

using namespace ProtoMol;
using namespace std;
using namespace ProtoMolAddon::Constant;

const string OutputToF::keyword("ToF");

OutputToF::OutputToF(): Output(-1)
{}

OutputToF::~OutputToF() {}

OutputToF::OutputToF(const string& filename):
  Output(1),
  reader(filename),
  detector(reader.GetValue<double>("detector.pos")),
  output_filename(reader.GetValue<string>("output_filename"))
{
}

Output* OutputToF::doMake(const std::vector<Value> &values) const
{
  return new OutputToF(values[0]);
}

void OutputToF::getParameters(std::vector<Parameter> &parameter) const
{
  parameter.push_back(Parameter(getId(), Value(output_filename, ConstraintValueType::NotEmpty())));
}

     
void OutputToF::doInitialize() {
  detector.Initialize(app);
}

void OutputToF::doRun(int step) {
  detector.UpdateRecord(app);
}

void OutputToF::doFinalize(int step) {
  ofstream f(output_filename);

  if (!f)
    std::cerr << "Fail to open file " << output_filename << std::endl;

  f << detector;
}

