#ifndef _REACTION_H_
#define _REACTION_H_

#include <random>
#include <vector>
#include <string>
#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/ProtoMolApp.h>

namespace ProtoMolAddon {
  namespace Reaction {

    template <class Mechanism>
    class Reaction {

    private:
      Mechanism m;
      unsigned int size;
      std::vector<Util::SIAtomProxy> atom_proxy;
      std::vector<typename Mechanism::state> atom_state;


    public:
      Reaction() : m(), size(0), atom_proxy(), atom_state() {}
      
      Reaction(const std::string &fname) : 
	m(fname), size(0), atom_proxy(), atom_state() {}

      void init(ProtoMol::ProtoMolApp *app) {
	size = app->positions.size();

	for (unsigned int i=0; i<size; i++) {
	  atom_proxy.push_back(Util::SIAtomProxy(app, i));
	  atom_state.push_back(m.assign_init_state(atom_proxy.back()));
	}
      }

      void react(double dt) {
	for (unsigned int i=0; i<size; i++)
	  m.react(atom_proxy[i], atom_state[i], dt);
      }
    };
  }
}

#endif
