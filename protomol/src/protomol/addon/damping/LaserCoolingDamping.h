#ifndef _LASER_COOLING_DAMPING_H_
#define _LASER_COOLING_DAMPING_H_

#include <protomol/type/Vector3D.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/util/SIAtomProxy.h>
#include <vector>
#include <string>
#include <iostream>

namespace ProtoMolAddon {
  namespace Damping {
    
    using namespace ProtoMol;
    
    class LaserCoolingDamping {

    private:
      struct Spec {
	
	struct Beam {
	  std::string label;
	  std::string ion_name;
	  double t_start;
	  double t_end;
	  Vector3D n;
	  double k;  // in cm^-1, w/o 2p
	  double s;   
	  double delta; // in Hz
	  double gamma; // in Hz

	  Beam() {}

	  Beam(const std::string &label, const std::string &ion_name, double t_start, double t_end, 
		const Vector3D &n_, double k, double s, double delta, double gamma) :
	    label(label), ion_name(ion_name), t_start(t_start), t_end(t_end), n(n_), k(k),
	    s(s), delta(delta), gamma(gamma) {
	    
	    n.normalize();
	  }
	  
	  Vector3D GetForce(const Vector3D &vel) const {
	    double d_eff = (delta - n * vel * k * 1e2)/gamma;
	    double pop = s / ( 1 + s + 4 * d_eff * d_eff) / 2;
	    
	    return n * (k * 1e2 * pop * Constant::HBAR * 2 * M_PI * gamma );
	  }

	  friend ostream& operator<< (ostream &os, Beam &entry) {
	    os << entry.ion_name << std::endl
	       << entry.t_start << std::endl
	       << entry.t_end << std::endl
	       << entry.n << std::endl
	       << entry.s << std::endl
	       << entry.delta << std::endl
	       << entry.gamma << std::endl;

	    return os;	    
	  }
	  
	  friend bool operator< (const Beam &e1, const Beam &e2) { return e1.ion_name < e2.ion_name; }
	};

	std::vector<Beam> beam_list;
	
	Spec() {}
	Spec(const std::string &fname);
      };
      
    private:
      Spec spec;

    public:
      LaserCoolingDamping(const LaserCoolingDamping::Spec &spec = LaserCoolingDamping::Spec());
      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;

      static std::string GetName() { return "LaserCoolingDampingForce"; }
      static std::string GetParameterName() { return "-laser-cooling-damping-spec"; }
    };

  }
}

#endif

