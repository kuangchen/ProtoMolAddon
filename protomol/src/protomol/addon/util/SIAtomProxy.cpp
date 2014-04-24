#include <protomol/addon/util/SIAtomProxy.h>

using namespace ProtoMolAddon::Util;

SIAtomProxy::SIAtomProxy() :
  name(NULL), mass(NULL), charge(NULL), pos(NULL), vel(NULL) {
}

SIAtomProxy::SIAtomProxy(ProtoMolApp *app, unsigned int i) :
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
  name(&(topo->atoms[i].name)),
  mass(&(topo->atoms[i].scaledMass)),
  charge(&(topo->atoms[i].scaledCharge)),
  pos(&(*positions)[i]),
  vel(&(*velocities)[i])
{  
}

ConstSIAtomProxy::ConstSIAtomProxy(const ProtoMolApp *app, unsigned int i) :
  name(&(app->topology->atoms[i].name)),
  mass(&(app->topology->atoms[i].scaledMass)),  
  charge(&(app->topology->atoms[i].scaledCharge)),
  pos(&(app->positions[i])),
  vel(&(app->velocities[i])) {
}

ConstSIAtomProxy::ConstSIAtomProxy(const GenericTopology *topo, 
				   const Vector3DBlock *positions, 
				   const Vector3DBlock *velocities, 
				   unsigned int i) :
  name(&(topo->atoms[i].name)),
  mass(&(topo->atoms[i].scaledMass)),
  charge(&(topo->atoms[i].scaledCharge)),
  pos(&(*positions)[i]),
  vel(&(*velocities)[i])
{  
}
