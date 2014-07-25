#include <protomol/addon/reaction/Ionization.h>
#include <cmath>
#include <protomol/type/Vector3D.h>

namespace ProtoMolAddon {
  namespace Reaction {

    namespace pt = boost::property_tree;
    namespace algorithm = boost::algorithm;

    Ionization::Spec::Spec(const std::string &fname):
      entry_list() {

      pt::ptree tree;
      pt::read_xml(fname, tree);

      r0 = tree.get<double>("ConfigRoot.IonizationSpec.r0");
      z0 = tree.get<double>("ConfigRoot.IonizationSpec.z0");

      for (auto &v: tree.get_child("ConfigRoot.IonizationSpec")) 
	if (v.first=="entry") 
	  entry_list.push_back(entry(v));
    }

    const std::string Ionization::keyword("Ionization");

    Ionization::Ionization(const std::string &fname) : 
      engine(rd()), dist(0, 1), spec(fname)
    {
    }
    
    Ionization::state Ionization::assign_init_state(Util::SIAtomProxy &ap) const {
      return find_entry(ap) == spec.entry_list.end() ?
	state::non_target : 
	state::neutral;
    }

    void Ionization::react(Util::SIAtomProxy &ap, Ionization::state &s, double dt) const {
      // if escaped or not a target, don't do anything
      //if (s!=state::non_target && s!=state::escaped)
      if (s==state::ionized) {
	// check if the ion has escaped
	ProtoMol::Vector3D p(ap.GetPosition());
	if (((p[0] * p[0] + p[0] * p[0]) > spec.r0 * spec.r0) || 
	    p[2] * p[2] > spec.z0 * spec.z0)
	  s = state::escaped;
      }

      else {
	spec_entry_iterator it = find_entry(ap);
	if (it != spec.entry_list.end())
	  if (exp(-it->rate * dt) < dist(engine)) {
	    std::cout << "atom " << ap.GetName() << " is ionized" << std::endl;
	    ap.SetIntegralCharge(it->charge);
	    s = state::ionized;
	  }
      }
    }

    
    Ionization::spec_entry_iterator Ionization::find_entry(Util::SIAtomProxy &ap) const { 
      if (ap.GetIntegralCharge()!=0)
	return spec.entry_list.end();

      auto it = std::find_if(spec.entry_list.begin(),
			     spec.entry_list.end(),
			     [&ap](const spec_entry &e) { return e.atom_name==ap.GetName(); });

      return it;
    }

  }
}
