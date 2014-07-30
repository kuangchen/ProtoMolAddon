#ifndef __ESCAPE_DETECTION_H
#define __ESCAPE_DETECTION_H

#include <string>
#include <protomol/addon/util/SIAtomProxy.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>

namespace ProtoMolAddon {
  namespace Util {
    class SIAtomProxy;
  }

  namespace Reaction {
    namespace pt = boost::property_tree;
    namespace algorithm = boost::algorithm;

    class EscapeDetection {
    private:

      struct Spec {
	double r0, z0;

	Spec(): r0(0), z0(0) {}
	Spec(const pt::ptree &t):
	  r0(t.get<double>("ConfigRoot.EscapeDetection.r0")),
	  z0(t.get<double>("ConfigRoot.EscapeDetection.z0")) 
	  {}

      };

    private:
      Spec spec;

    public:
      static const std::string keyword;
      enum state { not_escaped, escaped };
      
      EscapeDetection() {}
      EscapeDetection(const std::string &fname);

      state assign_init_state(Util::SIAtomProxy &ap) const;
      void react(Util::SIAtomProxy &ap, state &s, double dt) const;
    };
  }
}


#endif
