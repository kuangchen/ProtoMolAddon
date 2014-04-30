#ifndef _SIMPLE_DAMPING_H
#define _SIMPLE_DAMPING_H

#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/type/Vector3D.h>
#include <vector>
#include <string>
#include <iostream>

namespace ProtoMolAddon {
  namespace Damping {

    using namespace ProtoMol;
    
    class SimpleDamping  {

    private:
      struct Spec {

	struct Entry {
	  std::string label;
	  std::string ion_name;
	  double t_start;
	  double t_end;
	  double alpha;

	  Entry() {}

	  Entry(const std::string &label, const std::string &ion_name, 
		double t_start, double t_end, double alpha) : 
	    label(label), ion_name(ion_name), t_start(t_start), t_end(t_end), alpha(alpha) {};

	  friend bool operator< (const Entry &e1, const Entry &e2) { return e1.ion_name < e2.ion_name; }
	}; 

	std::vector<Entry> entry_list;

	Spec() : entry_list(0) {}
	Spec(const std::string &fname);
      };

    private:
      Spec spec;

    public:
      SimpleDamping();
      SimpleDamping(const Spec &spec);
      
      static std::string GetName() { return "SimpleDampingForce"; }
      static std::string GetParameterName() { return "-simple-damping-spec"; }
      Vector3D GetForce(const Util::ConstSIAtomProxy &atom, double now) const;
    };
  }
}

#endif
