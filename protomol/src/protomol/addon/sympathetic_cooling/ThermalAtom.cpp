#include <protomol/addon/sympathetic_cooling/ThermalAtom.h>
#include <protomol/base/PMConstants.h>
#include <cmath>

using namespace ProtoMolAddon::SympatheticCooling;
using namespace ProtoMol;
using namespace ProtoMol::Constant;

ThermalAtom::ThermalAtom(double m, double density, double temperature) : 
  mass(m), position(), velocity(), density(density), temperature(temperature),
  rd(), dice(0, sqrt(SI::BOLTZMANN*temperature/mass))
{

}  

void ThermalAtom::Resample() {
  velocity = Vector3D(dice(rd), dice(rd), dice(rd));
}

