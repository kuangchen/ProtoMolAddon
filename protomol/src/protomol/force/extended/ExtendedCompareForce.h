/* -*- c++ -*- */
#ifndef EXTENDEDCOMPAREFORCE_H
#define EXTENDEDCOMPAREFORCE_H

#include <protomol/force/CompareForce.h>
#include <protomol/force/extended/ExtendedForce.h>

namespace ProtoMol {
  //________________________________________ ExtendedCompareForce

  class ExtendedCompareForce : public CompareForce, public ExtendedForce {
    // This class contains the definition of one force

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ExtendedCompareForce(Force *actualForce, CompareForce *compareForce);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class ExtendedForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions,
                          const Vector3DBlock *velocities,
                          Vector3DBlock *forces,
                          ScalarStructure *energies);

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  const Vector3DBlock *velocities,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies);
  };

  //________________________________________ INLINES
}
#endif /* EXTENDEDCOMPAREFORCE_H */
