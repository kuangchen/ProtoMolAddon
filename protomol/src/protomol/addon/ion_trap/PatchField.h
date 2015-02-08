#ifndef _PATCH_FIELD_H
#define _PATCH_FIELD_H

#include <protomol/type/Vector3D.h>
#include <iostream>
#include <string>

namespace ProtoMolAddon {
  namespace Util {
    class ConstSIAtomProxy;
  }

  namespace IonTrap {
    using namespace ProtoMol;

    class PatchField {
    public:
      struct Spec {
	std::string label;
	double r0;
	double v;
	
	Spec() : label(""), r0(0), v(0) {}
	Spec(const std::string &fname);

	friend std::ostream& operator<< (std::ostream &os, const Spec &spec) {
	  os << "Patch Field Force" << std::endl 
	     << "label = " << spec.label << std::endl
	     << "r0 = " << spec.r0 << std::endl
	     << "v = " << spec.v << std::endl;

	  return os;
	}
      };

    private:
      Spec spec;
      
    public:
      PatchField(const Spec &spec = Spec());
      
      static std::string GetName() { return "PatchFieldForce"; }
      static std::string GetParameterName() { return "-pf-spec"; }
      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;

    };
  }
}

#endif

