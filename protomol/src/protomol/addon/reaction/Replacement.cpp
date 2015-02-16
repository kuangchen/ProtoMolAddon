#include <protomol/addon/reaction/Replacement.h>
#include <protomol/type/Vector3D.h>
#include <cmath>
#include <algorithm>
#include <boost/property_tree/xml_parser.hpp>

namespace ProtoMolAddon {
  namespace Reaction {

    namespace pt = boost::property_tree;

    Replacement::Spec::Spec() {}

    Replacement::Spec::Spec(const std::string &fname) {
      pt::ptree tree;
      pt::read_xml(fname, tree);
      
      for (auto &v: tree.get_child("ConfigRoot.ReplacementSpec"))
	if (v.first=="entry") 
	  (*this)[v.second.get<std::string>("from_name")]
	    = Replacement::rule(v.second.get<double>("rate"),
				v.second.get<std::string>("to_name"),
				v.second.get<double>("to_mass"),
				v.second.get<double>("to_energy"));

    }

    Replacement::Replacement() {}

    Replacement::Replacement(const Replacement& other):
      spec(other.spec),
      ap_array_ptr(other.ap_array_ptr ?
			 new Util::SIAtomProxyArray(*other.ap_array_ptr): NULL),
      engine(rd()),
      dist(0, 1)
    {
    }
    
    void Replacement::Initialize(ProtoMolApp *app) {
      ap_array_ptr.reset(new Util::SIAtomProxyArray(app));
    }

    void Replacement::Update(double now, double dt) {

      for (auto &ap: (*ap_array_ptr)) {
	const string &name = ap.GetName();
	
	if (spec.count(name) > 0) {
	  const rule &entry = spec[name];

	  if (exp(-entry.rate * dt) < dist(engine)) {
	    Vector3D v;
	    for (int j=0; j<3; j++) v[j] = dist(engine);

	    v.normalize();

	    double mag = sqrt(entry.to_energy * 2 / (entry.to_mass * Constant::ToSI::mass));
	    ap.SetVelocity(v * mag);	    
	    ap.SetName(entry.to_name);
	    ap.SetMass(entry.to_mass * Constant::ToSI::mass);
	  }
	}
      }
    }
  }
}
