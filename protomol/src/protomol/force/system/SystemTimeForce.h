/* -*- c++ -*- */
#ifndef SYSTEMTIMEFORCE_H
#define SYSTEMTIMEFORCE_H

#include <protomol/force/TimeForce.h>
#include <protomol/force/system/SystemForce.h>

namespace ProtoMol {
  //________________________________________ SystemTimeForce

  class SystemTimeForce : public TimeForce, public SystemForce {
    // This class contains the definition of one force

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    SystemTimeForce(Force *actualForce);
    virtual ~SystemTimeForce() {};
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
#endif /* SYSTEMTIMEFORCE_H */
