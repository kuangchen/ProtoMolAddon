#ifndef _EVAPORATION_H
#define _EVAPORATION_H

#include <string>
#include <memory>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/util/SIAtomProxyArray.h>

namespace ProtoMolAddon {
  namespace Util {
    class SIAtomProxyArray;
  }
  
  namespace Reaction {

    using namespace ProtoMol;
    
    class Evaporation {
    private:
      struct Spec {
	double r0;
	double z0;
	Spec() : r0(0), z0(0) {}
	Spec(const std::string &fname);
      };

      Spec spec;
      std::unique_ptr<Util::SIAtomProxyArray> ap_array_ptr;
      
    public:
      Evaporation() {}
      Evaporation(const Spec &spec) : spec(spec) {}
      Evaporation(const Evaporation &other);
      
    public:
      void Initialize(ProtoMolApp *app);
      void Update(double now, double dt);
      
      static const std::string GetName() { return "Evaporation"; }
      static const std::string GetParameterName() { return "-evaporation-spec"; }
    };
  }
}


#endif
