#include <protomol/module/IntegratorBaseModule.h>

#include <protomol/integrator/base/LangevinImpulseIntegrator.h>
#include <protomol/integrator/base/LangevinLeapfrogIntegrator.h>
#include <protomol/integrator/base/CGMinimizerIntegrator.h>
#include <protomol/integrator/base/NumericalDifferentiation.h>
#include <protomol/integrator/base/NumericallyDifferentiatedHessian.h>
#include <protomol/integrator/base/RMTIntegrator.h>

#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol;

void IntegratorBaseModule::init(ProtoMolApp *app) {
  app->integratorFactory.registerExemplar(new LangevinImpulseIntegrator());
  app->integratorFactory.registerExemplar(new LangevinLeapfrogIntegrator());
  app->integratorFactory.registerExemplar(new CGMinimizerIntegrator());
  app->integratorFactory.registerExemplar(new NumericalDifferentiation());
  app->integratorFactory.registerExemplar(new NumericallyDifferentiatedHessian());
  app->integratorFactory.registerExemplar(new RMTIntegrator());  

}
