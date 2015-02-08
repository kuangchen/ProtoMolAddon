#include <protomol/addon/damping/SimpleDamping.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>

namespace ProtoMolAddon {
  namespace Damping {

    namespace pt = boost::property_tree;
    namespace algorithm = boost::algorithm;

    SimpleDamping::Spec::Spec(const std::string &fname) 
      : entry_list() 
    {
      typedef SimpleDamping::Spec::Entry entry;
      pt::ptree tree;
      pt::read_xml(fname, tree);
 
      for (auto &v : tree.get_child("ConfigRoot.SimpleDampingSpec")) 
	if (v.first == "entry") {
	  std::string target_atom(v.second.get<std::string>("target_atom"));
	  algorithm::trim(target_atom);

	  entry_list.push_back(entry(v.second.get<std::string>("label"),
				     target_atom,
				     v.second.get<double>("t_start"),
				     v.second.get<double>("t_end"),
				     v.second.get<double>("alpha")));
  
	}

      std::sort(entry_list.begin(), entry_list.end());
    }

    SimpleDamping::SimpleDamping() {}

    SimpleDamping::SimpleDamping(const SimpleDamping::Spec &spec) :
      spec(spec) {}


    Vector3D SimpleDamping::GetForce(const Util::ConstSIAtomProxy &atom, double now) const {
      typedef SimpleDamping::Spec::Entry spec_entry;

      if (spec.entry_list.empty()) return Vector3D();

      auto iter = std::lower_bound(spec.entry_list.begin(),
				   spec.entry_list.end(),
				   atom.GetName(),
				   [](const spec_entry &e, const std::string &s) {
				     return e.target_atom < s;
				   });

      if (iter == spec.entry_list.end() ) 
	return Vector3D();

      else 
	return (iter->t_start > now || iter->t_end < now) ? Vector3D() : atom.GetVelocity() * (-iter->alpha);
    }


    std::ostream& operator<< (std::ostream &os, const SimpleDamping& d) {
      std::string header("target_atom\tt_start\tt_end\talpha\tlabel\n");
      
      os << header;
      for (auto &e: d.spec.entry_list)
	os << e.target_atom << "\t" 
	   << e.t_start << "\t" 
	   << e.t_end << "\t" 
	   << e.alpha << "\t" 
	   << e.label << std::endl;
      
      return os;
    }
    
  }
}
