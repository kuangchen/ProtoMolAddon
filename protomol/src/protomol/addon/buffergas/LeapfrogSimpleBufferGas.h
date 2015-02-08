#ifndef _LEAPFROG_SIMPLE_BUFFERGAS_H_
#define _LEAPFROG_SIMPLE_BUFFERGAS_H_

#include <memory>
#include <protomol/type/Vector3D.h>
#include <random>
#include <string>

namespace ProtoMol {
  class ProtoMolApp;
}

namespace ProtoMolAddon {
  namespace Util {
    class SIAtomProxy;
    class SIAtomProxyArray;
  }
  
  namespace BufferGas {

    using namespace ProtoMol;
    
    class LeapfrogSimpleBufferGas {

      struct NeutralAtom {
      private:
	double m;
	std::string name;
	double alpha;
	double T;
	double rho;
	std::string target_atom_name;

	std::random_device rd;
	std::default_random_engine engine;
	std::normal_distribution<double> vel_dist;
	std::exponential_distribution<double> interval;
	std::uniform_real_distribution<double> uniform_dist;
	double sigma;
	Vector3D v;

      private:
	double GetCollisionInterval(Util::SIAtomProxy &ap);
	
      public:
	NeutralAtom() {};
	NeutralAtom(const std::string &fname);
	void Collide(Util::SIAtomProxy &ap, double dt);
      };

    private:
      NeutralAtom neutral;
      std::unique_ptr<Util::SIAtomProxyArray> ap_array_ptr;

    public:
      LeapfrogSimpleBufferGas();
      LeapfrogSimpleBufferGas(const std::string &fname);

      static std::string GetName() { return "LeapFrogSimpleBufferGas"; }
      static std::string GetParameterName() { return "-simple-buffergas-spec"; }

      void Initialize(ProtoMolApp *app);
      void Update(double now, double dt);
    };

  }
}


#endif
