#include <protomol/addon/util/ConstSIAtomProxy.h>
#include <protomol/ProtoMolApp.h>

namespace ProtoMolAddon {
  namespace Util {
    
    using namespace ProtoMol;
    
    ConstSIAtomProxy::ConstSIAtomProxy():
      id(-1), name(NULL), mass(NULL), charge(NULL), pos(NULL), vel(NULL) {}

    ConstSIAtomProxy::ConstSIAtomProxy(const ConstSIAtomProxy &other) {
      id = other.id;
      name = other.name;
      mass = other.mass;
      charge = other.charge;
      pos = other.pos;
      vel = other.vel;
    }
    
    ConstSIAtomProxy::ConstSIAtomProxy(const ProtoMolApp *app, unsigned int i):
      id(i),
      name(&(app->topology->atoms[i].name)),
      mass(&(app->topology->atoms[i].scaledMass)),  
      charge(&(app->topology->atoms[i].scaledCharge)),
      pos(&(app->positions[i])),
      vel(&(app->velocities[i]))
    {}

    ConstSIAtomProxy::ConstSIAtomProxy(const GenericTopology *topo, 
			     const Vector3DBlock *positions, 
			     const Vector3DBlock *velocities, 
			     const unsigned int i) :
      id(i),
      name(&(topo->atoms[i].name)),
      mass(&(topo->atoms[i].scaledMass)),
      charge(&(topo->atoms[i].scaledCharge)),
      pos(&(*positions)[i]),
      vel(&(*velocities)[i])
    {}
  }
}
