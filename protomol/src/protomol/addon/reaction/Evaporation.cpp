#include <protomol/addon/reaction/Evaporation.h>
#include <boost/property_tree/xml_parser.hpp>

namespace ProtoMolAddon {
  namespace Reaction {

    namespace pt = boost::property_tree;

    Evaporation::Spec::Spec(const std::string &fname) {
      pt::ptree tree;
      pt::read_xml(fname, tree);

      r0 = tree.get<double>("ConfigRoot.Evaporation.r0");
      z0 = tree.get<double>("ConfigRoot.Evaporation.z0");
    }
    
    void Evaporation::Initialize(ProtoMolApp *app) {
      ap_array_ptr.reset(new Util::SIAtomProxyArray(app));
    }

    void Evaporation::Update(double now, double dt) {
      for (auto &ap: *ap_array_ptr) {
	const Vector3D &pos = ap.GetPosition();

	if ( (ap.GetIntegerCharge()!=0) &&
	     ((pos[0]*pos[0] + pos[1]*pos[1]) > spec.r0 * spec.r0 ||
	      (fabs(pos[2]) > spec.z0)) )
	    ap.SetIntegerCharge(0);
      }
    }
  }
}
    


