#include <protomol/force/extended/ExtendedCompareForce.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/topology/GenericTopology.h>

using namespace ProtoMol;
//____ ExtendedCompareForce

ExtendedCompareForce::ExtendedCompareForce(Force *actualForce,
                                           CompareForce *compareForce) :
  CompareForce(actualForce, compareForce) {}

void ExtendedCompareForce::evaluate(const GenericTopology *topo,
                                    const Vector3DBlock *positions,
                                    const Vector3DBlock *velocities,
                                    Vector3DBlock *forces,
                                    ScalarStructure *energies) {
  preprocess(positions->size());
  (dynamic_cast<ExtendedForce *>(myActualForce))->evaluate(topo,
    positions,
    velocities,
    myForces,
    myEnergies);
  postprocess(topo, forces, energies);
}

void ExtendedCompareForce::parallelEvaluate(const GenericTopology *topo,
                                            const Vector3DBlock *positions,
                                            const Vector3DBlock *velocities,
                                            Vector3DBlock *forces,
                                            ScalarStructure *energies) {
  preprocess(positions->size());
  (dynamic_cast<ExtendedForce *>(myActualForce))->parallelEvaluate(topo,
    positions,
    velocities,
    myForces,
    myEnergies);
  postprocess(topo, forces, energies);
}

