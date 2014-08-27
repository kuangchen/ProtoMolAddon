#ifndef _STRAY_FIELD_H
#define _STRAY_FIELD_H

#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/type/Vector3D.h>
#include <iosfwd>

namespace ProtoMolAddon {
  namespace Util {
    class ConstSIAtomProxy;
  }

  namespace StrayField {

    using namespace ProtoMol;

    class StrayField {
    public:
      struct Spec {
	Vector3D field;
	Spec() {}
	Spec(const std::string &fname);
      };

    private:
      Spec spec;

    public:
      StrayField() {}
      StrayField(const Spec &spec);

      static std::string GetName() { return "StrayFieldForce"; }
      static std::string GetParameterName() { return "-stray-field-spec"; }
      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;
      friend std::ostream& operator<< (std::ostream &os, const StrayField &f);
    };

    
  }
}
    
 

#endif

