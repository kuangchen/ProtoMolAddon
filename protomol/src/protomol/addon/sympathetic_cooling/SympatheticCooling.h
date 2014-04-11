#ifndef __SYMPATHETIC_COOLING_H
#define __SYMPATHETIC_COOLING_H

#include <protomol/addon/util/AtomProxy.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/ProtoMolApp.h>
#include <vector>
#include <memory>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace ProtoMol;
    using namespace std;

    template <class Atom, class Collision>
    class SympatheticCooling {
    private:
      Atom atom;
      Collision collision;
      random_device rd;
      uniform_real_distribution<double> dice;

    public:
      SympatheticCooling(ProtoMolApp *app);
      //void Collide(double dt) const;
    };


    template <class Atom, class Collision>
    SympatheticCooling<Atom, Collision>::SympatheticCooling(ProtoMolApp *app): 
      atom(), collision(),
      rd(), dice(0, 1),
    {

    }


    template <class Atom, class Collision>
    void SympatheticCooling<Atom, Collision>::Collide(double dt) const {

      for (auto &ion : ion_proxy) {
	double rate = collsion.GetReactionRateConstant(ion) * atom.density;

	if (1-exp(-rate * dt) > dice(rd)) {
	  atom.Resample();
	  
	  double b1 = ion->mi / (atom.mass + ion->mi);
	  double b2 = 1-b1;
	  
	  Vector3D v_i(ion->GetVelocity()), v_n(atom.velocity);
	  Vector3D v_rel(v_i, v_n);
	  Vector3D v_com(b1*v_i + b2*v_n);

	  Vector3D v_rel_after = collision.Rotate(v_rel);
	  ion.SetVelocity(v_com + v_rel_after * b2);
	}
      }

    }
  }
}

#endif
