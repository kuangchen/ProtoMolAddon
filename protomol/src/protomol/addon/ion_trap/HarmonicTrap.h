#ifndef _HARMONIC_TRAP_H_
#define _HARMONIC_TRAP_H_

#include <string>
#include <protomol/type/Vector3D.h>

namespace ProtoMolAddon {
  namespace Util {
    class ConstSIAtomProxy;
  }

  namespace IonTrap {
    using namespace ProtoMol;

    class HarmonicTrap {
    private:
      struct Spec {
	Vector3D omega;
	Spec() {}
	Spec(const std::string &fname);
      };

    private:
      Spec spec;

    public:
      HarmonicTrap(const Spec &spec = Spec());

      static std::string GetName() { return "HarmonicTrapForce"; }
      static std::string GetParameterName() { return "-ht-spec"; }
      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;
    };
    
  }
}

#endif
