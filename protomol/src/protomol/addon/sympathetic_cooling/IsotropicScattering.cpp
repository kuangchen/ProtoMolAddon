#include <protomol/addon/sympathetic_cooling/IsotropicScattering.h>
#include <cmath>

using namespace ProtoMolAddon::SympatheticCooling;

IsotropicScattering::IsotropicScattering(const initializer &init) :
  rd(), z_dice(-1, 1), phi_dice(0, M_PI*2)
{

}

Vector3D IsotropicScattering::Rotate(double mu, const Vector3D &v) const {
  double mag = v.norm();

  double z = z_dice(rd);
  double phi = phi_dice(rd);
  double r = sqrt(1 - z*z);
  double x = cos(phi) * r;
  double y = sin(phi) * r;

  return Vector3D(x, y, z) * mag;
}

double IsotropicScattering::GetReactionRateConstant(const ThermalAtom& atom, const SIAtomProxy &ion) const {
  double mu = atom.GetMass() * ion.GetMass() / (atom.GetMass() + ion.GetMass());

  return 2*M_PI * sqrt(atom.GetC4()/mu);
}
