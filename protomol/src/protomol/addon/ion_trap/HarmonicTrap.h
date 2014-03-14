#ifndef _HARMONIC_TRAP_H_
#define _HARMONIC_TRAP_H_

#include <iosfwd>
#include <array>
#include <protomol/type/Vector3D.h>

using ProtoMol::Vector3D;
using std::istream;
using std::ostream;
using std::array;

namespace ProtoMolAddon {
  namespace IonTrap {

    class HarmonicTrap {
    private:
      array<double, 3> freq;

    public:
      HarmonicTrap(double omega_x=0, double omega_y=0, double omega_z=0);
      Vector3D GetForce(const Vector3D &pos, double mass) const;

      friend istream& operator>> (istream &is, HarmonicTrap &t);
      friend ostream& operator<< (ostream &os, const HarmonicTrap &t);
    };
    
  }
}

#endif
