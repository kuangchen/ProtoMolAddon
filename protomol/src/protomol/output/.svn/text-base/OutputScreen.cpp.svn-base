#include <protomol/output/OutputScreen.h>
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

const string OutputScreen::keyword("Screen");


OutputScreen::OutputScreen() : unit("fs"), factor(1.0) {}


OutputScreen::OutputScreen(int freq) : Output(freq), unit("fs"), factor(1) {}


void OutputScreen::doInitialize() {
  Real step = app->integrator->getTimestep() *
    max(1, getOutputFreq());

  if (step >= 1e13) {
    unit = "s";
    factor = 1e-15;

  } else if (step >= 1e10) {
    unit = "ms";
    factor = 1e-12;

  } else if (step >= 1e7) {
    unit = "us";
    factor = 1e-9;

  } else if (step >= 1e4) {
    unit = "ns";
    factor = 1e-6;

  } else if (step >= 1e1) {
    unit = "ps";
    factor = 1e-3;
  }
}


void OutputScreen::doRun(int step) {
  report << plain << "Step : ";
  report.setf(ios::right);
  report.width(7);
  report << step << ", Time : ";
  report.width(10);
  report.setf(ios::showpoint | ios::fixed);
  report.precision(3);
  report << app->outputCache.getTime() * factor << " [" << unit << "], TE : ";
  report.precision(4);
  report.width(12);
  report << app->outputCache.getTotalEnergy() << " [kcal/mol]";
  report << ", T : ";
  report.precision(4);
  report.width(10);
  report << app->outputCache.getTemperature() << " [K]";
  report << ", V : ";
  report.precision(2);
  report.width(12);
  report << app->outputCache.getVolume() << " [AA^3]" << endr;
  report.reset();
}


Output *OutputScreen::doMake(const vector<Value> &values) const {
  return new OutputScreen(values[1]);
}


bool OutputScreen::isIdDefined(const Configuration *config) const {
  return config->valid("outputFreq") && !config->empty(getId()) &&
    (!config->valid(getId()) || ((*config)[getId()] == true));
}


void OutputScreen::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back(Parameter(getId(), Value(true), true));
  Output::getParameters(parameter);
}
