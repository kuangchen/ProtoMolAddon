/* -*- c++ -*- */
#ifndef SYSTEMCOMPAREFORCE_H
#define SYSTEMCOMPAREFORCE_H

#include <protomol/force/CompareForce.h>
#include <protomol/force/system/SystemForce.h>
#include <protomol/force/hessian/ReducedHessAngle.h>

namespace ProtoMol {
  //________________________________________ SystemCompareForce

  class SystemCompareForce : public CompareForce, public SystemForce {
    // This class contains the definition of one force

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    SystemCompareForce(Force *actualForce, CompareForce *compareForce);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    using SystemForce::evaluate; // Avoid compiler warning/error
    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions,
                          Vector3DBlock *forces,
                          ScalarStructure *energies);

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies);
  };

  //________________________________________ INLINES
}
#endif /* SYSTEMCOMPAREFORCE_H */
