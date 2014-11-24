\#ifndef __GENERIC_COLLISION_H
#define __GENERIC_COLLISION_H

#include <protomol/type/Vector3D.h>
#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/addon/sympathetic_cooling/ThermalAtom.h>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace ProtoMolAddon::Util;
    using namespace ProtoMol;
    using namespace std;

    class GenericScattering {
    public:
      GenericScattering() {};
      virtual Vector3D Rotate(double mu, const Vector3D &v) const = 0;
      virtual double GetReactionRateConstant(const ThermalAtom &atom, const SIAtomProxy &ion) const = 0;
    };

  }
}

#endif
