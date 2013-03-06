/* -*- c++ -*- */
#ifndef MOLLYFORCE_H
#define MOLLYFORCE_H

#include <protomol/force/Force.h>
#include <vector>

namespace ProtoMol {
  class ForceGroup;
  class GenericTopology;
  class Vector3DBlock;
  class ScalarStructure;
  class ReducedHessAngle;

  //________________________________________ MollyForce

  class MollyForce : public Force {
    // This class contains the definition of one Molly force (a force that
    // works only on positions)

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    MollyForce() : Force() {}
    virtual ~MollyForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class MollyForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    virtual void evaluate(const GenericTopology *,
                          const Vector3DBlock *,
                          std::vector<ReducedHessAngle> *) = 0;
    // Evaluate this force using the information given.


    virtual void parallelEvaluate(const GenericTopology *,
                                  const Vector3DBlock *,
                                  std::vector<ReducedHessAngle> *);
    // Evaluate this force using the information given.
    // This is a default methods for forces which do not have their
    // own parallel version, but are called from a parallel environment.
    // A parallel environment is when the total forces and energies
    // is the sum over all nodes, where in a sequential environment all
    // nodes have the total forces and energies locally.
    // Since parallelEvaluate() is called we are in a parallel environment,
    // so only one node will evaluate the force. This gives us a parallelization
    // for free (for bonded forces).
    //
    // To the implementor:
    // Overload this methods to implement a parallel version. Split up the
    // interactions and distribute them among the nodes. This can by done
    // either by ranges or sequences. Overload numberOfBlocks() such that the
    // number matches the next() calls inside parallelEvaluate().

    virtual void addToForceGroup(ForceGroup *forceGroup);
  };

  //________________________________________ INLINES

  inline void MollyForce::parallelEvaluate(
    const GenericTopology *topo,
    const Vector3DBlock *positions,
    std::vector<ReducedHessAngle> *
    angleFilter) {
    evaluate(topo, positions, angleFilter);
  }
}
#endif /* MOLLYFORCE_H */



