#include <protomol/module/IntegratorOpenMMModule.h>

#ifdef HAVE_OPENMM
#include <protomol/integrator/openMM/OpenMMIntegrator.h>
#endif

#ifdef HAVE_OPENMM_LTMD
#include <protomol/integrator/openMM/NormalModeOpenMM.h>
#endif

#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol;

void IntegratorOpenMMModule::init(ProtoMolApp *app) {
#ifdef HAVE_OPENMM
  app->integratorFactory.registerExemplar(new OpenMMIntegrator());
#endif

#ifdef HAVE_OPENMM_LTMD
  app->integratorFactory.registerExemplar(new NormalModeOpenMM());
#endif
}
