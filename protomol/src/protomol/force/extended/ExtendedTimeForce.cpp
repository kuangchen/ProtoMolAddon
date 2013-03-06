#include <protomol/force/extended/ExtendedTimeForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/topology/GenericTopology.h>

using namespace ProtoMol;
//____ ExtendedTimeForce

ExtendedTimeForce::ExtendedTimeForce(Force *actualForce) :
  TimeForce(
    actualForce) {}

void ExtendedTimeForce::evaluate(const GenericTopology *topo,
                                 const Vector3DBlock *positions,
                                 const Vector3DBlock *velocities,
                                 Vector3DBlock *forces,
                                 ScalarStructure *energies) {
  preprocess(positions->size());
  (dynamic_cast<ExtendedForce *>(myActualForce))->evaluate(topo,
    positions,
    velocities,
    forces,
    energies);
  postprocess(topo, forces, energies);
}

void ExtendedTimeForce::parallelEvaluate(const GenericTopology *topo,
                                         const Vector3DBlock *positions,
                                         const Vector3DBlock *velocities,
                                         Vector3DBlock *forces,
                                         ScalarStructure *energies) {
  preprocess(positions->size());
  (dynamic_cast<ExtendedForce *>(myActualForce))->parallelEvaluate(topo,
    positions,
    velocities,
    forces,
    energies);
  postprocess(topo, forces, energies);
}

