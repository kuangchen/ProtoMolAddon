#include <protomol/addon/damping/SimpleDamping.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <stdexcept>

namespace pt = boost::property_tree;
using namespace ProtoMolAddon::Damping;
namespace algorithm = boost::algorithm;

SimpleDamping::Spec::Spec(const std::string &fname) 
  : entry_list() 
{
  typedef SimpleDamping::Spec::Entry entry;
  pt::ptree tree;
  pt::read_xml(fname, tree);
 
  for (auto &v : tree.get_child("ConfigRoot.SimpleDampingSpec")) 
    if (v.first == "Entry") {

      // Need to strip the whitespace in atom_name
      std::string atom_name(v.second.get<std::string>("atom_name"));
      algorithm::trim(atom_name);

      entry_list.push_back(entry(v.second.get<std::string>("label"),
				 atom_name,
				 v.second.get<double>("t_start"),
				 v.second.get<double>("t_end"),
				 v.second.get<double>("alpha")));
  
    }

  std::sort(entry_list.begin(), entry_list.end());

  std::cout << (*this);
}

SimpleDamping::SimpleDamping() {}

SimpleDamping::SimpleDamping(const SimpleDamping::Spec &spec) :
  spec(spec) {}


Vector3D SimpleDamping::GetForce(const Util::ConstSIAtomProxy &atom, double now) const {
  typedef SimpleDamping::Spec::Entry spec_entry;

  if (spec.entry_list.empty()) 
    return Vector3D();

  auto iter = std::lower_bound(spec.entry_list.begin(),
  			       spec.entry_list.end(),
  			       atom.GetName(),
  			       [](const spec_entry &e, const std::string &s) {
  				 return e.atom_name < s;
  			       });

  if (iter == spec.entry_list.end() ) 
    return Vector3D();

  else 
    return (iter->t_start > now || iter->t_end < now) ? Vector3D() : atom.GetVelocity() * (-iter->alpha);
  
}

namespace ProtoMolAddon {
  namespace Damping {

    std::ostream& operator<< (std::ostream &os, const SimpleDamping& d) {
      std::string header("atom_name\tt_start\tt_end\talpha\tlabel\n");
      
      os << header;
      for (auto &e: d.spec.entry_list)
	os << e.atom_name << "\t" 
	   << e.t_start << "\t" 
	   << e.t_end << "\t" 
	   << e.alpha << "\t" 
	   << e.label << std::endl;
      
      return os;
    }
  }
}
