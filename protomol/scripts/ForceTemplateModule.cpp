#include "TemplateModule.h"

#include <protomol/base/ProtoMolApp.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/TopologyModule.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/VacuumBoundaryConditions.h>

using namespace std;
using namespace ProtoMol;

void TemplateModule::registerForces(ProtoMolApp *app) {
  ForceFactory &f = app->forceFactory;
  string boundConds = app->config[InputBoundaryConditions::keyword];

  if (equalNocase(boundConds, PeriodicBoundaryConditions::keyword)) {
  } else if (equalNocase(boundConds, VacuumBoundaryConditions::keyword)) {
  }
}
