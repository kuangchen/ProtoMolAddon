#include <protomol/addon/tof/CEM.h>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>

namespace ProtoMolAddon {
  namespace ToF {

    namespace pt = boost::property_tree;
    namespace algorithm = boost::algorithm;
    
    CEM::Spec::Spec(const std::string &conf_fname) {
      pt::ptree tree;
      pt::read_xml(conf_fname, tree);

      x_cem = tree.get<double>("ConfigRoot.CEM.x_cem");
      r_cem = tree.get<double>("ConfigRoot.CEM.r_cem");
      fname = tree.get<std::string>("ConfigRoot.CEM.fname");
      algorithm::trim(fname);
    }

    CEM::CEM(const CEM::Spec &spec) :
      spec(spec) {
    }

    void CEM::Initialize(const ProtoMolApp *app) {
      ap_array_ptr.reset(new Util::ConstSIAtomProxyArray(app));
    }

    void CEM::Update(double now) {
      for (auto ap: *ap_array_ptr) {
	if (hit_entry_map.count(ap) == 0) {
	  const Vector3D &pos = ap.GetPosition();
	  if (pos[0] > spec.x_cem &&
	      ((pos[1]*pos[1] + pos[2]*pos[2]) < spec.r_cem * spec.r_cem))
	    
	    hit_entry_map[ap] = HitEntry(now, pos, ap.GetVelocity());
	}
      }
    }

    void CEM::Finalize() {
      std::ofstream os(spec.fname);
      
      for (auto &kv: hit_entry_map) {
	os << kv.first.GetName() << "\t"
	   << kv.second.t << "\t"
	   << kv.second.pos << "\t"
	   << kv.second.vel << std::endl;
      }
    }
  }
}

