#ifndef __SYMPATHETIC_COOLING_H_
#define __SYMPATHETIC_COOLING_H_

#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/ProtoMolApp.h>
#include <vector>

namespace ProtoMolAddon {
  namespace SympatheticCooling {

    using namespace ProtoMolAddon::Util;
    using namespace ProtoMol;
    using namespace std;

    template <class Atom, class Collision>
    class SympatheticCooling {
    private:
      mutable Atom atom;
      Collision collision;
      vector<SIAtomProxy> ion_proxy_list;

      random_device rd;
      mutable uniform_real_distribution<double> dice;

    public:
      SympatheticCooling();
      SympatheticCooling(const typename Atom::Spec &atom_spec, 
			 const typename Collision::initializer &collision_init);
      void Initialize(ProtoMolApp *app);
      void Run(double dt);
    };

    template <class Atom, class Collision>
    SympatheticCooling<Atom, Collision>::SympatheticCooling() {}

    template <class Atom, class Collision>
    SympatheticCooling<Atom, Collision>::SympatheticCooling(const typename Atom::Spec &atom_spec,
							    const typename Collision::initializer &collision_init): 
      atom(atom_spec), collision(collision_init),
      ion_proxy_list(),
      rd(), dice(0, 1)
    {
    }


    template <class Atom, class Collision>
    void SympatheticCooling<Atom, Collision>::Initialize(ProtoMolApp *app) {
      for (unsigned int i=0; i<app->positions.size(); i++) 
	ion_proxy_list.push_back(SIAtomProxy(app, i));
    }


    template <class Atom, class Collision>
    void SympatheticCooling<Atom, Collision>::Run(double dt) {

      for (auto &ion : ion_proxy_list) {
	double rate = collision.GetReactionRateConstant(atom, ion) * atom.GetDensity();
	
	if (1-exp(-rate * dt) > dice(rd)) {
	  atom.Resample();
	  
	  double b1 = ion.GetMass() / (atom.GetMass() + ion.GetMass());
	  double b2 = 1-b1;
	  
	  Vector3D v_i(ion.GetVelocity()), v_n(atom.GetVelocity());
	  Vector3D v_rel(v_i - v_n);
	  Vector3D v_com(v_i*b1 + v_n*b2);
	  
	  double mu = ion.GetMass() * atom.GetMass() / (ion.GetMass()+atom.GetMass());
	  Vector3D v_rel_after = collision.Rotate(mu, v_rel);
	  ion.SetVelocity(v_com + v_rel_after * b2);
	}
      }

    }
  }
}

#endif
