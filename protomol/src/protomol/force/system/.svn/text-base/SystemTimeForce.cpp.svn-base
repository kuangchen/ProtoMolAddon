#include <protomol/force/system/SystemTimeForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/topology/GenericTopology.h>

using namespace std;
using namespace ProtoMol;
//____ SystemTimeForce

SystemTimeForce::SystemTimeForce(Force *actualForce) :
  TimeForce(actualForce)
{}

void SystemTimeForce::evaluate(const GenericTopology *topo,
                               const Vector3DBlock *positions,
                               Vector3DBlock *forces,
                               ScalarStructure *energies) {
  preprocess(positions->size());
  (dynamic_cast<SystemForce *>(myActualForce))->evaluate(topo,
    positions,
    forces,
    energies);
  postprocess(topo, forces, energies);
}

void SystemTimeForce::parallelEvaluate(const GenericTopology *topo,
                                       const Vector3DBlock *positions,
                                       Vector3DBlock *forces,
                                       ScalarStructure *energies) {
  preprocess(positions->size());
  (dynamic_cast<SystemForce *>(myActualForce))->parallelEvaluate(topo,
    positions,
    forces,
    energies);
  postprocess(topo, forces, energies);
}

