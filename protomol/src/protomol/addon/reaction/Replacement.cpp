#include <protomol/addon/Constants.h>
#include <protomol/addon/reaction/Replacement.h>
#include <cmath>
#include <algorithm>
#include <protomol/type/Vector3D.h>
#include <protomol/addon/util/SIAtomProxy.h>


namespace ProtoMolAddon {
  namespace Reaction {

    namespace pt = boost::property_tree;
    namespace algorithm = boost::algorithm;


    Replacement::Spec::Spec(const std::string &fname):
      entry_list() {

      pt::ptree tree;
      pt::read_xml(fname, tree);

      for (auto &v: tree.get_child("ConfigRoot.ReplacementSpec"))
	if (v.first=="entry") entry_list.push_back(spec_entry(v));
    }

    const std::string Replacement::keyword("Replacement");

    Replacement::Replacement(const std::string &fname) : 
      engine(rd()), dist(0, 1), spec(fname)
    {}
    
    Replacement::state Replacement::assign_init_state(Util::SIAtomProxy &ap) const {
      return (find_entry(ap) == spec.entry_list.end()) ? 
	state::non_target :
	state::not_replaced;
    }

    void Replacement::react(Util::SIAtomProxy &ap, Replacement::state &s, double dt) const {
      // if replaced  a target, don't do anything
      //if (s!=state::non_target && s!=state::escaped)
      if (s==state::not_replaced) {
	spec_entry_iterator it = find_entry(ap);
	if (it != spec.entry_list.end()) 
	    if (exp(-it->rate * dt) < dist(engine)) {
	      std::cout << "atom " << ap.GetName() << " is replaced to " << it->to_name << std::endl;

	      ProtoMol::Vector3D v;
	      for (int j=0; j<3; j++) 
		v[j] = dist(engine);

	      v.normalize();

	      double mag = sqrt(it->to_energy * 2 / (it->to_mass * Constant::ToSI::mass));
	      ap.SetVelocity(v * mag);
	      ap.SetName(it->to_name);
	      ap.SetMass(it->to_mass * Constant::ToSI::mass);
	      s = state::replaced;
	    }
      }
    }

    Replacement::spec_entry_iterator Replacement::find_entry(Util::SIAtomProxy &ap) const { 
      auto it = std::find_if(spec.entry_list.begin(),
			     spec.entry_list.end(),
			     [&ap](const spec_entry &e) { return e.from_name==ap.GetName(); });

      return it;
    }

  }
}
