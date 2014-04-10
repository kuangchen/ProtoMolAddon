#include <protomol/addon/sympathetic_cooling/ProtoMolIonProxy.h>
#include <protomol/ProtoMolApp.h>

using namespace ProtoMolAddon;
using namespace ProtoMol::Constant;
using namespace ProtoMolAddon::SympatheticCooling;

ProtoMolIonProxy::ProtoMolIonProxy(ProtoMol::ProtoMolApp *app, unsigned int i):
  mass(app->topology->atoms[i].scaledMass),
  charge(app->topology->atoms[i].scaledCharge),
  position(app->positions[i]),
  velocity(app->velocities[i])
{}




