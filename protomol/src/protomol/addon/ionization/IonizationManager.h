#ifndef __IONIZATION_MANAGER_H_
#define __IONIZATION_MANAGER_H_

#include <random>
#include <vector>
#include <string>
#include <protomol/addon/util/SIAtomProxy.h>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Ionization {
    
    class IonizationManager {
      
      struct Spec {
	struct entry {
	  std::string atom_name;
	  double ionization_rate;
	  int charge;
	  
	  entry(const std::string &name, double rate, int charge): 
	    atom_name(name), ionization_rate(rate), charge(charge) {}
	};

	std::vector<entry> entry_list;

      public:
	Spec(): entry_list() {}	
	Spec(const std::string &fname);

	friend std::ostream& operator<< (std::ostream &os, const Spec &spec) {
	  for (auto &e: spec.entry_list)
	    os << e.atom_name << "\t" << e.ionization_rate << "\t" << e.charge << std::endl;
	  return os;
	}

      };

      
    private:
      Spec spec;
      std::vector<Util::SIAtomProxy> atom_proxy_list;
      std::map<unsigned int, unsigned int> proxy_map;
      std::random_device rd;
      std::default_random_engine engine;
      std::uniform_real_distribution<> dist;
      
    private:
      void Ionize(double dt);

    public:
      IonizationManager() {}
      IonizationManager(const std::string &fname) : spec(fname), proxy_map(), engine(rd()), dist(0, 1) {}
      void Init(ProtoMol::ProtoMolApp *app);
    };
  }
}

#endif
