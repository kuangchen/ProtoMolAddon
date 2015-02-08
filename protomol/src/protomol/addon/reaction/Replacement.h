#ifndef _REPLACEMENT_H
#define _REPLACEMENT_H

#include <random>
#include <map>
#include <string>
#include <memory>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/util/SIAtomProxyArray.h>

namespace ProtoMolAddon {
  namespace Util {
    class SIAtomProxyArray;
  }
  
  namespace Reaction {
    using namespace ProtoMol;
    
    class Replacement {
    private:
      struct rule {
	double rate;
	std::string to_name;
	double to_mass;
	double to_energy;
	
	rule(): rate(0), to_name(), to_mass(0), to_energy(0) {}
	rule(double rate, const std::string& to_name,
	     double to_mass, double to_energy):
	  rate(rate),
	  to_name(to_name),
	  to_mass(to_mass),
	  to_energy(to_energy)
	{} 
	  
      };

      struct Spec: public std::map<std::string, rule> {
      public:
	Spec();
	Spec(const std::string &fname);
      };
	
    private:
      Spec spec;
      std::unique_ptr<Util::SIAtomProxyArray> ap_array_ptr;

      std::random_device rd;
      mutable std::default_random_engine engine;
      mutable std::uniform_real_distribution<> dist;

    public:
      Replacement();
      Replacement(const Spec& spec): spec(spec), engine(rd()), dist(0, 1) {}

      void Initialize(ProtoMolApp *app);
      void Update(double now, double dt);

      static const std::string GetName() { return "Replacement"; }
      static const std::string GetParameterName() { return "-replacement-spec"; }
    };
  }
}


#endif
