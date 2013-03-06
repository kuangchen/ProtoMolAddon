#include "IntegratorSCEModule.h"

#include "../integrator/LangevinFlowCoupledIntegrator.h"

#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol;

void IntegratorSCEModule::init(ProtoMolApp *app) {
  app->integratorFactory.registerExemplar(new LangevinFlowCoupledIntegrator());

}
