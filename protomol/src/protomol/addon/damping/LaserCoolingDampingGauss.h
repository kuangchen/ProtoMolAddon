#ifndef _LASER_COOLING_DAMPING_GAUSS_H_
#define _LASER_COOLING_DAMPING_GAUSS_H_

#include <protomol/type/Vector3D.h>
#include <protomol/addon/Constants.h>
#include <protomol/addon/util/ConstSIAtomProxy.h>
#include <protomol/addon/util/SIAtomProxy.h>
#include <vector>
#include <string>
#include <iostream>

namespace ProtoMolAddon {
  namespace Damping {
    
    using namespace ProtoMol;
    
    class LaserCoolingDampingGauss {

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
	  double waist; // in m

	  Beam() {}

	  Beam(const std::string &label, const std::string &ion_name, double t_start, double t_end, 
		const Vector3D &n_, double k_, double s_, double delta_, double gamma_, double waist_) :
	    label(label), ion_name(ion_name), t_start(t_start), t_end(t_end), n(n_), k(k_*1e2 * 2 * M_PI),
	    s(s_), delta(delta_ * 2 * M_PI), gamma(gamma_ * 2 * M_PI), waist(waist_) {
	    
	    n.normalize();
	    
	  }
	  
	  Vector3D GetForce(const Vector3D &vel, const Vector3D &pos) const {
	    double p = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );
            double l_n = sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] );
	    double n_0 = n[0]/l_n; double n_1 = n[1]/l_n; double n_2 = n[2]/l_n;
	    double dot = pos[0] * n_0 + pos[1] * n_1 + pos[2] * n_2;
	    double r = sqrt( p * p - dot * dot );
	    double d_eff = (delta - n * vel * k)/gamma;
	    double s_exp = s * exp ( -2 * ( r * r ) / ( waist * waist ));
	    double pop = s_exp / ( 1 + s_exp + 4 * d_eff * d_eff) / 2;
	    return n * (k * Constant::HBAR * pop * gamma );
	  }

	  friend ostream& operator<< (ostream &os, Beam &entry) {
	    os << entry.ion_name << std::endl
	       << entry.t_start << std::endl
	       << entry.t_end << std::endl
	       << entry.n << std::endl
	       << entry.s << std::endl
	       << entry.k << std::endl
	       << entry.delta << std::endl
	       << entry.gamma << std::endl
               << entry.waist << std::endl;

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
      LaserCoolingDampingGauss(const LaserCoolingDampingGauss::Spec &spec = LaserCoolingDampingGauss::Spec());
      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;

      static std::string GetName() { return "LaserCoolingDampingGaussForce"; }
      static std::string GetParameterName() { return "-laser-cooling-damping-gauss-spec"; }
    };

  }
}

#endif

