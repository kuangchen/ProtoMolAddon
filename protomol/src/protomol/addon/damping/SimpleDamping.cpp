#include <protomol/addon/damping/SimpleDamping.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <algorithm>
#include <stdexcept>

namespace pt = boost::property_tree;
using namespace ProtoMolAddon::Damping;

SimpleDamping::Spec::Spec(const std::string &fname) 
  : entry_list() 
{
  typedef SimpleDamping::Spec::Entry entry;
  pt::ptree tree;
  pt::read_xml(fname, tree);
 
  for (auto &v : tree.get_child("SimpleDampingSpec")) 
    if (v.first == "Entry") 
      entry_list.push_back(entry(v.second.get<std::string>("name"),
				 v.second.get<double>("t_start"),
				 v.second.get<double>("t_end"),
				 v.second.get<double>("alpha")));
  

  std::sort(entry_list.begin(), entry_list.end());
}

SimpleDamping::SimpleDamping() {}

SimpleDamping::SimpleDamping(const SimpleDamping::Spec &spec) :
  spec(spec) {}


Vector3D SimpleDamping::GetForce(const Util::ConstSIAtomProxy &atom, double now) const {
  typedef SimpleDamping::Spec::Entry spec_entry;

  if (spec.entry_list.empty())
    return Vector3D();

  std::vector<spec_entry>::const_iterator iter = std::lower_bound(spec.entry_list.begin(),
								  spec.entry_list.end(),
								  atom.GetName(),
								  [](const spec_entry &e, const std::string &s) {
								    return e.name < s;
								  });
  if (iter == spec.entry_list.end() )
    return Vector3D();
  
  else 
    return (iter->t_start > now || iter->t_end < now) ? Vector3D() : atom.GetVelocity() * (-iter->alpha);
  
}
