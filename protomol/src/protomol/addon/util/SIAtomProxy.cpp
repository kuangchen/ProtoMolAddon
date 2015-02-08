#include <protomol/addon/util/SIAtomProxy.h>
#include <protomol/ProtoMolApp.h>

namespace ProtoMolAddon {
  namespace Util {

    SIAtomProxy::SIAtomProxy() :
      id(-1), name(NULL), mass(NULL), charge(NULL), pos(NULL), vel(NULL) {
    }

    SIAtomProxy::SIAtomProxy(const SIAtomProxy &other) {
      id = other.id;
      name = other.name;
      mass = other.mass;
      charge = other.charge;
      pos = other.pos;
      vel = other.vel;
    }
    
    
    SIAtomProxy::SIAtomProxy(ProtoMolApp *app, unsigned int i) :
      id(i),
      name(&(app->topology->atoms[i].name)),
      mass(&(app->topology->atoms[i].scaledMass)),  
      charge(&(app->topology->atoms[i].scaledCharge)),
      pos(&(app->positions[i])),
      vel(&(app->velocities[i])) {

    }

    SIAtomProxy::SIAtomProxy(GenericTopology *topo, 
			     Vector3DBlock *positions, 
			     Vector3DBlock *velocities, 
			     unsigned int i) :
      id(i),
      name(&(topo->atoms[i].name)),
      mass(&(topo->atoms[i].scaledMass)),
      charge(&(topo->atoms[i].scaledCharge)),
      pos(&(*positions)[i]),
      vel(&(*velocities)[i])
    {  
    }

  }
}
