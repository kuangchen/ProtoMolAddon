#ifndef _RADIAL_QUENCHING_H
#define _RADIAL_QUENCHING_H

#include <protomol/type/Vector3D.h>
#include <iostream>
#include <string>

namespace ProtoMolAddon {
  namespace Util {
    class ConstSIAtomProxy;
  }

  namespace IonTrap {
    using namespace ProtoMol;

    class RadialQuenching {
    public:
      struct Spec {
	std::string label;
	double r0;
	double v;
	
	Spec() : label(""), r0(0), v(0) {}
	Spec(const std::string &fname);

	friend std::ostream& operator<< (std::ostream &os, const Spec &spec) {
	  os << "Radial Quenching Force" << std::endl 
	     << "label = " << spec.label << std::endl
	     << "r0 = " << spec.r0 << std::endl
	     << "v = " << spec.v << std::endl;

	  return os;
	}
      };

    private:
      Spec spec;
      
    public:
      RadialQuenching(const Spec &spec = Spec());
      
      static std::string GetName() { return "RadialQuenchingForce"; }
      static std::string GetParameterName() { return "-radial-quenching-spec"; }
      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;

    };
  }
}

#endif

