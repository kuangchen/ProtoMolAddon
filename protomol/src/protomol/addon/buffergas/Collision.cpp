#include <protomol/ProtoMolApp.h>
#include <protomol/addon/buffergas/Collision.h>
#include <protomol/addon/util/SIAtomProxy.h>

namespace ProtoMolAddon {
  namespace BufferGas {

    Collision::Collision() : atom_proxy() {};

    void Collision::Initialize(ProtoMol::ProtoMolApp *app) {
      for (unsigned int i=0; i<app->positions.size(); i++)
	atom_proxy.push_back(Util::SIAtomProxy(app, i));
    }

    void Collision::Collide(double dt) {
      for (auto &ap: atom_proxy) 
	CollideEach(ap, dt);
    }

  }
}
