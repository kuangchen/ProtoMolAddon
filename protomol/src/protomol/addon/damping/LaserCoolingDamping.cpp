#include <protomol/addon/damping/LaserCoolingDamping.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <algorithm>
#include <stdexcept>

namespace pt = boost::property_tree;
using namespace ProtoMolAddon::Damping;


LaserCoolingDamping::Spec::Spec(const std::string &fname) {
  typedef LaserCoolingDamping::Spec::Entry entry;
  pt::ptree tree;
  pt::read_xml(fname, tree);
 
  for (auto &v : tree.get_child("LaserCoolingDampingSpec")) 
    if (v.first == "Entry") 
      entry_list.push_back(entry(v.second.get<std::string>("ion_name"),
				 v.second.get<double>("t_start"),
				 v.second.get<double>("t_end"),
				 v.second.get<Vector3D>("n"),
				 v.second.get<double>("k"),
				 v.second.get<double>("s"),
				 v.second.get<double>("delta"),
				 v.second.get<double>("gamma")));
  
  std::sort(entry_list.begin(), entry_list.end());
}

LaserCoolingDamping::LaserCoolingDamping(const LaserCoolingDamping::Spec &spec) :
  spec(spec) {
}

Vector3D LaserCoolingDamping::GetForce(const Util::ConstSIAtomProxy &atom, double now) const {
  typedef LaserCoolingDamping::Spec::Entry entry;

  if (spec.entry_list.empty())
    return Vector3D();

  std::vector<entry>::const_iterator iter = std::lower_bound(spec.entry_list.begin(),
							     spec.entry_list.end(),
							     atom.GetName(),
							     [](const entry &e, const std::string &s) {
							       return e.ion_name < s;
							     });

  if (iter == spec.entry_list.end() )
    return Vector3D();
  
  else 
    return (iter->t_start > now || iter->t_end < now) ? Vector3D() : iter->GetForce(atom.GetVelocity());
}

