#ifndef _CEM_H
#define _CEM_H

#include <protomol/type/Vector3DBlock.h>
#include <protomol/addon/util/ConstSIAtomProxyArray.h>
#include <memory>
#include <map>
#include <string>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {

  namespace Util {
    class ConstSIAtomProxy;
    class ConstSIAtomProxyArray;
  }
  
  namespace ToF {
    
    using namespace ProtoMol;

    class CEM {
    public:
      typedef enum ion_status { flying = 0, hit = 1 } ion_status;
      
      struct Spec {
	double x_cem;
	double r_cem;
	std::string fname;

	Spec() {}
	Spec(const std::string &fname);
      };

      struct HitEntry {
	double t;
	Vector3D pos;
	Vector3D vel;
	HitEntry(double t, Vector3D pos, Vector3D vel) :
	  t(t), pos(pos), vel(vel) {}

	HitEntry() {}
      };
      
    private:
      Spec spec;
      std::map<Util::ConstSIAtomProxy, HitEntry> hit_entry_map;
      std::unique_ptr<Util::ConstSIAtomProxyArray> const_ap_array_ptr;
      
    public:
      CEM() {}
      CEM(const CEM &other);
      CEM(const Spec &spec);
      
      void Initialize(const ProtoMolApp *app);
      void Update(double now);
      void Finalize();

      static std::string GetName() { return "CEM"; }
      static std::string GetParameterName() { return "-cem-spec"; }
    };

  }
}


#endif
