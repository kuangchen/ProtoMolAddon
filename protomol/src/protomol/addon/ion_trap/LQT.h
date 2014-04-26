#ifndef _LQT_H
#define _LQT_H

#include <protomol/type/Vector3D.h>
#include <string>

namespace ProtoMolAddon {
  namespace Util {
    class ConstSIAtomProxy;
  }

  namespace IonTrap {
    using ProtoMol::Vector3D;

    class LQT {

      struct LQTSpec {
	double r0;
	double z0;
	double omega;
	double v_rf;
	double v_ec;
	double kappa;

	LQTSpec() {};
	LQTSpec(const std::string &fname);
      };

    private:
      LQTSpec spec;
      double cache_a, cache_b, cache_c;
      
      //MathieuFunc::mathieu_param GetMathieuParam(const Util::ConstSIAtomProxy &atom) const;
      
    public:
      LQT();
      LQT(const LQTSpec &spec);

      //array<double, 3> GetSecularFrequency(const Util::ConstSIAtomProxy &atom) const;

      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;

    };
  }
}


#endif
