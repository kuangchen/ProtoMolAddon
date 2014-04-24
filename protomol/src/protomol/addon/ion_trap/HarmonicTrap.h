#ifndef _HARMONIC_TRAP_H_
#define _HARMONIC_TRAP_H_

#include <iosfwd>
#include <string>
#include <array>
#include <protomol/type/Vector3D.h>
#include <protomol/addon/util/SIAtomProxy.h>

using ProtoMol::Vector3D;
using namespace std;

namespace ProtoMolAddon {
  namespace IonTrap {

    class HarmonicTrap {

    private:
      struct HarmonicTrapSpec {
	double omega[3];
	HarmonicTrapSpec() {}
	HarmonicTrapSpec(const string &fname);
      };

    private:
      HarmonicTrapSpec spec;

    public:
      HarmonicTrap();
      HarmonicTrap(const HarmonicTrapSpec &spec);

      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;
    };
    
  }
}

#endif
