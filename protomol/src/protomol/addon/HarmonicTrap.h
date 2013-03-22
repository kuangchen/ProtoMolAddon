#ifndef _HARMONIC_TRAP_H
#define _HARMONIC_TRAP_H

#include <protomol/type/Vector3D.h>
#include <protomol/base/PMConstants.h>
#include <vector>

using namespace std;
using namespace ProtoMol;

namespace ProtoMolAddon {
  
  class HarmonicTrap {
  private:
    vector<double> freq;

  public:
    HarmonicTrap();
    HarmonicTrap(const string &def);
    void GetForce(const Vector3D &pos, double mass, Vector3D& force);
  };

}

#endif
