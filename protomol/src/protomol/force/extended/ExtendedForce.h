/* -*- c++ -*- */
#ifndef EXTENDEDFORCE_H
#define EXTENDEDFORCE_H

#include <protomol/force/Force.h>

namespace ProtoMol {
  class CompareForce;
  class TimeForce;
  class ForceGroup;
  class GenericTopology;
  class Vector3DBlock;
  class CompareForce;
  class TimeForce;

  //________________________________________ ExtendedForce

  class ExtendedForce : virtual public Force {
    // This class contains the definition of one extended force (a force 
    // that works on unprocessed positions and velocities)

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ExtendedForce() : Force() {}
    virtual ~ExtendedForce() {}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions,
                          const Vector3DBlock *velocities,
                          Vector3DBlock *forces,
                          ScalarStructure *energies) = 0;
    // Evaluate this force using the information given.

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  const Vector3DBlock *velocities,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies);

    virtual CompareForce *makeCompareForce(Force *actualForce,
                                           CompareForce *compareForce) const;
    virtual TimeForce *makeTimeForce(Force *actualForce) const;

    virtual void addToForceGroup(ForceGroup *forceGroup);
  };
  //________________________________________ INLINES

  inline void ExtendedForce::parallelEvaluate(const GenericTopology *topo,
                                              const Vector3DBlock *positions,
                                              const Vector3DBlock *velocities,
                                              Vector3DBlock *forces,
                                              ScalarStructure *energies) {
    evaluate(topo, positions, velocities, forces, energies);
  }
}
#endif /* EXTENDEDFORCE_H */
