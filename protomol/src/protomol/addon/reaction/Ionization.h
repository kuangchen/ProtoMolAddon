#ifndef __IONIZAION_H_
#define __IONIZAION_H_

#include <random>
#include <vector>
#include <string>
#include <protomol/addon/util/SIAtomProxy.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>


namespace ProtoMolAddon {
  namespace Util {
    class SIAtomProxy;
  }

  namespace Reaction {
    namespace pt = boost::property_tree;
    namespace algorithm = boost::algorithm;

    class Ionization {
    private:

      struct Spec {

	struct entry {
	  double rate;
	  std::string atom_name;
	  int charge;

	  entry(double rate, const std::string &atom_name, int charge):
	    rate(rate), atom_name(atom_name), charge(charge) {}

	  entry(const pt::ptree::value_type &v) :
	    rate(v.second.get<double>("rate")),
	    atom_name(v.second.get<std::string>("atom_name")),
	    charge(v.second.get<int>("to_energy")) 
	    {
	      algorithm::trim(atom_name);
	    }
	};

	std::vector<entry> entry_list;
	double r0, z0;

	Spec(): entry_list(), r0(0), z0(0) {}
	Spec(const std::string &fname);

      };

    private:
      typedef Spec::entry spec_entry;
      typedef std::vector<Spec::entry>::const_iterator spec_entry_iterator;

      std::random_device rd;
      mutable std::default_random_engine engine;
      mutable std::uniform_real_distribution<> dist;
      Spec spec;

    private:
      spec_entry_iterator find_entry(Util::SIAtomProxy &ap) const;

    public:
      static const std::string keyword;
      enum state { non_target, neutral, ionized, escaped };
      
      Ionization(): engine(rd()), dist(0, 1) {}
      Ionization(const std::string &fname);

      state assign_init_state(Util::SIAtomProxy &ap) const;
      void react(Util::SIAtomProxy &ap, state &s, double dt) const;
    };
  }
}


#endif
