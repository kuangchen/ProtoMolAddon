#include <protomol/module/LeapfrogModule.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/integrator/leapfrog/LeapfrogIntegrator.h>
#include <protomol/integrator/leapfrog/LeapfrogTruncatedShadow.h>
#include <protomol/integrator/leapfrog/DMDLeapfrogIntegrator.h>
#include <protomol/integrator/leapfrog/PLeapfrogIntegrator.h>
#include <protomol/integrator/leapfrog/NoseNVTLeapfrogIntegrator.h>
#include <protomol/integrator/leapfrog/LeapfrogDataAcquisition.h>
#include <protomol/integrator/leapfrog/GPU.h>

using namespace ProtoMol;

void LeapfrogModule::init(ProtoMolApp *app) {
  app->integratorFactory.registerExemplar(new LeapfrogIntegrator());
  app->integratorFactory.registerExemplar(new LeapfrogTruncatedShadow());
  app->integratorFactory.registerExemplar(new DMDLeapfrogIntegrator());
  app->integratorFactory.registerExemplar(new PLeapfrogIntegrator());
  app->integratorFactory.registerExemplar(new NoseNVTLeapfrogIntegrator());
  app->integratorFactory.registerExemplar(new LeapfrogDataAcquisition());
  app->integratorFactory.registerExemplar(new GPU());
}
