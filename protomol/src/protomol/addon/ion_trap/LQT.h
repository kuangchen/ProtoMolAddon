#ifndef _LQT_H
#define _LQT_H

#include <iostream>
#include <stdexcept>
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
	LQTSpec(double r0, double z0, double omega, double v_rf, double v_ec, double kappa):
	  r0(r0), z0(z0), omega(omega), v_rf(v_rf), v_ec(v_ec), kappa(kappa) {
	  
	  if (r0<0 || z0<0 || omega<0 || v_rf<0 || v_ec<0 || kappa<0)
	    throw std::runtime_error("Invalid input for LQT Spec");

	}

	LQTSpec(const std::string &fname);

	friend std::ostream& operator<< (std::ostream &os, const LQTSpec &spec) {
	  os << "r0 = " << spec.r0 << std::endl
	     << "z0 = " << spec.z0 << std::endl
	     << "omega = " << spec.omega << std::endl
	     << "v_rf = " << spec.v_rf << std::endl
	     << "v_ec = " << spec.v_ec << std::endl
	     << "kappa = " << spec.kappa << std::endl;
	  return os;
	}
	
      };

    private:
      LQTSpec spec;
      double cache_a, cache_b, cache_c;
      
      //MathieuFunc::mathieu_param GetMathieuParam(const Util::ConstSIAtomProxy &atom) const;
      
    public:
      LQT();
      LQT(const LQTSpec &spec);

      //array<double, 3> GetSecularFrequency(const Util::ConstSIAtomProxy &atom) const;

      static std::string GetName() { return "LQTForce"; }
      static std::string GetParameterName() { return "-lqt-spec"; }
      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;

    };
  }
}


#endif
