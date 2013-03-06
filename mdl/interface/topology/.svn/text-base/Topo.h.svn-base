#ifndef TOPO_H
#define TOPO_H

#include <protomol/topology/TopologyUtilities.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
void remLin(Vector3DBlock *velocities,
                                const GenericTopology *topo) {
	removeLinearMomentum(velocities,topo);
}


void remAng(const Vector3DBlock* positions, Vector3DBlock *velocities,
                                const GenericTopology *topo) {
	removeAngularMomentum(positions, velocities,topo);
	}
	}

#endif
