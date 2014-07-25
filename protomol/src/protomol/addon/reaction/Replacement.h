#ifndef __REPLACEMENT_H_
#define __REPLACEMENT_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <random>
#include <string>
#include <vector>

namespace ProtoMolAddon {
  namespace Util {
    class SIAtomProxy;
  }

  namespace Reaction {
    namespace pt = boost::property_tree;
    namespace algorithm = boost::algorithm;

    class Replacement {
    private:
      struct Spec {
	
	struct entry {
	  double rate;
	  std::string from_name;
	  std::string to_name;
	  double to_mass;
	  double to_energy;

	  entry(double rate, const std::string &from_name, const std::string &to_name, double to_mass, double to_energy): 
	    rate(rate), from_name(from_name), to_name(to_name), to_mass(to_mass), to_energy(to_energy) {}
	
	  entry(const pt::ptree::value_type &v):
	    rate(v.second.get<double>("rate")),
	    from_name(v.second.get<std::string>("from_name")),
	    to_name(v.second.get<std::string>("to_name")),
	    to_mass(v.second.get<double>("to_mass")),
	    to_energy(v.second.get<double>("to_energy")) 
	    {
	      
	      algorithm::trim(from_name);
	      algorithm::trim(to_name);
	    }
	};
	public:
	  std::vector<entry> entry_list;
	  
	public:
	  Spec(): entry_list() {}
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
      enum state { non_target, not_replaced, replaced };
      
      Replacement(): engine(rd()), dist(0, 1) {}
      Replacement(const std::string &fname);

      state assign_init_state(Util::SIAtomProxy &ap) const;
      void react(Util::SIAtomProxy &ap, state &s, double dt) const;

    };
  }
}

#endif
