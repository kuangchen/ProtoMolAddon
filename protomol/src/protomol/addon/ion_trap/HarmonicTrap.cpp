#include <protomol/addon/ion_trap/HarmonicTrap.h>
#include <protomol/type/Vector3D.h>
#include <iostream>

using std::istream;
using std::ostream;
using namespace ProtoMolAddon::IonTrap;

HarmonicTrap::HarmonicTrap(double omega_x, double omega_y, double omega_z) {
  freq[0] = omega_x;
  freq[1] = omega_y;
  freq[2] = omega_z;
}

Vector3D HarmonicTrap::GetForce(const Vector3D &pos, double mass) const {
  Vector3D force;
  for (int i=0; i<3; i++)
    force[i] = -mass * freq[i] * freq[i] * pos[i];

  return force;
}

namespace ProtoMolAddon {
  namespace IonTrap {
    istream& operator>> (istream &is, HarmonicTrap &t) {
      is >> t.freq[0] >> t.freq[1] >> t.freq[2];
      return is;
    }

    ostream& operator<< (ostream &os, const HarmonicTrap &t) {
      os << t.freq[0] << t.freq[1] << t.freq[2];
      return os;
    }
  }
}
