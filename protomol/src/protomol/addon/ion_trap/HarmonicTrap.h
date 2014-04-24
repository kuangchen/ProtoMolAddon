#ifndef _HARMONIC_TRAP_H_
#define _HARMONIC_TRAP_H_

#include <string>
#include <array>
#include <protomol/type/Vector3D.h>

namespace ProtoMolAddon {
  namespace Util {
    class ConstSIAtomProxy;
  }

  namespace IonTrap {
    
    using ProtoMol::Vector3D;

    class HarmonicTrap {

    private:
      struct HarmonicTrapSpec {
	std::array<double, 3> omega;
	HarmonicTrapSpec() {}
	HarmonicTrapSpec(const std::string &fname);
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
