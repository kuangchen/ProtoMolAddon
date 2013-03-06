#include "OutputScreenUnit.h"
#include <protomol/output/OutputCache.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/config/Configuration.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/module/MainModule.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/OutputModule.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputScreenUnit
const string OutputScreenUnit::keyword("Screen");

OutputScreenUnit::OutputScreenUnit() :
  Output(), myUnit("us"), myFactor(1.0) {}

OutputScreenUnit::OutputScreenUnit(int freq) :
  Output(freq), myUnit("us"), myFactor(1.0) {}

void OutputScreenUnit::doInitialize() {
  Real step = app->integrator->getTimestep() *
    max(1, getOutputFreq());

  if (step >= 1e6) {
    myUnit = "s";
    myFactor = 1e-6;
  }else if (step >= 1e3) {
    myUnit = "ms";
    myFactor = 1e-3;
  }
}

void OutputScreenUnit::doRun(int step) {
  report << plain << "Step : ";
  report.setf(ios::right);
  report.width(7);
  report << step << ", Time : ";
  report.width(10);
  report.setf(ios::showpoint | ios::fixed);
  report.precision(3);
  report << app->outputCache.time() * myFactor << " [" << myUnit << "], TE : ";
  report.precision(4);
  report.width(12);
  report << app->outputCache.totalEnergy() << " [fJ]";
  report << ", T : ";
  report.precision(4);
  report.width(10);
  report << app->outputCache.temperature() << " [K]";
  report << ", V : ";
  report.precision(2);
  report.width(12);
  report << app->outputCache.volume() << " [um^3]" << endr;
  report.reset();
}

Output *OutputScreenUnit::doMake(const vector<Value> &values) const {
  return new OutputScreenUnit(values[1]);
}

bool OutputScreenUnit::isIdDefined(const Configuration *config) const {
  return config->valid("outputFreq") && !config->empty(getId()) &&
    (!config->valid(getId()) || ((*config)[getId()] == true));
}

void OutputScreenUnit::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back(Parameter(getId(), Value(true), true));
  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(getOutputFreq(), ConstraintValueType::Positive())));
}

bool OutputScreenUnit::adjustWithDefaultParameters(vector<Value> &values,
                                               const Configuration *config)
const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  return checkParameters(values);
}
