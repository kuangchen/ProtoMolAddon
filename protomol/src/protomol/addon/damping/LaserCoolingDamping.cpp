#include <protomol/addon/damping/LaserCoolingDamping.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <algorithm>
#include <stdexcept>

namespace pt = boost::property_tree;
using namespace ProtoMolAddon::Damping;

LaserCoolingDamping::Spec::Spec(const std::string &fname) {
  pt::ptree tree;
  pt::read_xml(fname, tree);
 
  for (auto &v : tree.get_child("ConfigRoot.LaserCoolingDampingSpec")) 
    if (v.first == "beam") {
      Beam b(v.second.get<std::string>("label"),
	     v.second.get<std::string>("ion_name"),
	     v.second.get<double>("t_start"),
	     v.second.get<double>("t_end"),
	     v.second.get<Vector3D>("n"),
	     v.second.get<double>("k"),
	     v.second.get<double>("s"),
	     v.second.get<double>("delta"),
	     v.second.get<double>("gamma"));

      beam_list.push_back(b);
    }

  // sort the beam based on ion_name
  // to make searching fast
  std::sort(beam_list.begin(), beam_list.end());
}

LaserCoolingDamping::LaserCoolingDamping(const LaserCoolingDamping::Spec &spec) :
  spec(spec) {
}

Vector3D LaserCoolingDamping::GetForce(const Util::ConstSIAtomProxy &atom, double now) const {
  typedef LaserCoolingDamping::Spec::Beam beam;

  if (spec.beam_list.empty())
    return Vector3D();
  
  auto lb = std::lower_bound(spec.beam_list.cbegin(), 
   			     spec.beam_list.cend(), 
   			     atom.GetName(), 
   			     [](const beam &e, const std::string &s) { return e.ion_name < s; });

  auto ub = std::upper_bound(lb, spec.beam_list.cend(),
   			     atom.GetName(),
   			     [](const std::string &s, const beam &e) { return s < e.ion_name; });

  // return zero force, if no beam is associated with the atom
  if (lb-ub == 0) return Vector3D();
  else {
    // otherwise add all the force associated with the atom
    Vector3D f;
	
    std::for_each(lb, ub, 
		  [&f, &now, &atom](const beam &e) {
		    //cout << e.label << "\n";
		    //for (double v=-10; v<10; v+=1) 
		    //  cout << e.GetForce(Vector3D(0, 0, v)) << "\n";

		    if (e.t_start < now && e.t_end > now) 
		      f += e.GetForce(atom.GetVelocity());
		  });
    return f;
  }
}

