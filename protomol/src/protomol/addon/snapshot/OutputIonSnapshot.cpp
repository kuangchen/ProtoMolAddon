#include <sstream>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/snapshot/OutputIonSnapshot.h>
#include <string>

using namespace ProtoMol;
using namespace ProtoMolAddon::Snapshot;
using namespace std;
using namespace ProtoMolAddon::Constant;

const string OutputIonSnapshot::keyword("IonSnapshot");

OutputIonSnapshot::OutputIonSnapshot(): Output(-1) {}

OutputIonSnapshot::~OutputIonSnapshot() {}

OutputIonSnapshot::OutputIonSnapshot(const string& fname)
  : Output(1),
    fname(fname),
    manager(fname)
{
}

Output* OutputIonSnapshot::doMake(const std::vector<Value> &values) const
{
  return new OutputIonSnapshot(values[0]);
}

void OutputIonSnapshot::getParameters(std::vector<Parameter> &parameter) const
{
  parameter.push_back(Parameter(getId(), Value(fname, ConstraintValueType::NotEmpty())));
}

     
void OutputIonSnapshot::doInitialize() {   
}

void OutputIonSnapshot::doRun(int step) {
  double now = app->topology->time * TIME_CONV;
  manager.Run(now, app);
}

void OutputIonSnapshot::doFinalize(int step) {
}

