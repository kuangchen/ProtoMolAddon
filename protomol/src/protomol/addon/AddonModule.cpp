#include <protomol/addon/AddonModule.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/LeapfrogBufferGasIntegrator.h>
#include <protomol/addon/LeapfrogBufferGasIntegrator2.h>
#include <protomol/addon/reaction/ReactionIntegrator.h>
#include <protomol/addon/reaction/Replacement.h>
#include <protomol/addon/reaction/Ionization.h>

namespace ProtoMolAddon {

  void AddonModule::init(ProtoMol::ProtoMolApp *app) {
    app->integratorFactory.registerExemplar(new LeapfrogBufferGasIntegrator());
    app->integratorFactory.registerExemplar(new LeapfrogBufferGasIntegrator2());
    app->integratorFactory.registerExemplar(new Reaction::ReactionIntegrator<Reaction::Replacement>);
    app->integratorFactory.registerExemplar(new Reaction::ReactionIntegrator<Reaction::Ionization>);
  }
}
