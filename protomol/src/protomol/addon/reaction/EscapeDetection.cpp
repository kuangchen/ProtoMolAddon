#include <protomol/addon/reaction/EscapeDetection.h>
#include <cmath>

namespace ProtoMolAddon {
  namespace Util {
    class SIAtomProxy;
  }

  namespace Reaction {

    const std::string EscapeDetection::keyword("EscapeDetection");

    EscapeDetection::EscapeDetection(const std::string &fname) {
      pt::ptree tree;
      pt::read_xml(fname, tree);

      spec = Spec(tree);
    }

    EscapeDetection::state EscapeDetection::assign_init_state(Util::SIAtomProxy &ap) const {
      return state::not_escaped;
    }

    void EscapeDetection::react(Util::SIAtomProxy &ap, EscapeDetection::state &s, double dt) const {
      ProtoMol::Vector3D p(ap.GetPosition());
      if ( (s==state::not_escaped) &&  
	   ((p[0]*p[0]+p[1]*p[1] > fabs(spec.r0)) || (fabs(p[2]) > fabs(spec.z0))) ) {
	
	ap.SetIntegralCharge(0);
	s = state::escaped;

	std::cout << p << " another one bites the dust\n";
      }
    
    }
  }
}


