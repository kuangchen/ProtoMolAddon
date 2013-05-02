#include <protomol/addon/ProtoMolIon.h>
#include <protomol/addon/Constants.h>

using namespace ProtoMolAddon;
using namespace ProtoMol::Constant;
using namespace ProtoMolAddon::Constant;

ProtoMolIon::ProtoMolIon(const ProtoMolApp *app, unsigned int i):
  mass(app->topology->atoms[i].scaledMass * SI::AMU),
  charge(app->topology->atoms[i].scaledCharge / SQRTCOULOMBCONSTANT * SI::ELECTRON_CHARGE),
  position(app->positions[i] * POSITION_CONV),
  velocity(app->velocities[i] * VELOCITY_CONV),
  atomId(i)
{
}

void ProtoMolIon::UpdateProtoMol(ProtoMolApp *app) {
  app->positions[atomId] = position / POSITION_CONV;
  app->velocities[atomId] = velocity / VELOCITY_CONV;
}


