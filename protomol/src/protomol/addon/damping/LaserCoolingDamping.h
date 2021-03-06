#ifndef _LASER_COOLING_DAMPING_H_
#define _LASER_COOLING_DAMPING_H_

#include <protomol/type/Vector3D.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/util/ConstSIAtomProxy.h>
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
	  std::string target_atom_name;
	  double t_start;
	  double t_end;
	  Vector3D n;
	  double k;  // in cm^-1, w/o 2p
	  double s;   
	  double delta; // in Hz
	  double gamma; // in Hz

	  Beam() {}

	  Beam(const std::string &label_, const std::string &target_atom_name, double t_start, double t_end, 
		const Vector3D &n_, double k_, double s_, double delta_, double gamma_) :
	    label(label_), target_atom_name(target_atom_name), t_start(t_start), t_end(t_end), n(n_), k(k_*1e2 * 2 * M_PI),
	    s(s_), delta(delta_ * 2 * M_PI), gamma(gamma_ * 2 * M_PI) {
	    
	    n.normalize();
	    
	  }
	  
	  Vector3D GetForce(const Vector3D &vel) const {
	    double d_eff = (delta - n * vel * k)/gamma;
	    double pop = s / ( 1 + s + 4 * d_eff * d_eff) / 2;
	    return n * (k * Constant::HBAR * pop * gamma );
	  }

	  friend ostream& operator<< (ostream &os, Beam &entry) {
	    os << entry.target_atom_name << std::endl
	       << entry.t_start << std::endl
	       << entry.t_end << std::endl
	       << entry.n << std::endl
	       << entry.s << std::endl
	       << entry.k << std::endl
	       << entry.delta << std::endl
	       << entry.gamma << std::endl;

	    return os;	    
	  }
	  
	  friend bool operator< (const Beam &e1, const Beam &e2) { return e1.target_atom_name < e2.target_atom_name; }
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

