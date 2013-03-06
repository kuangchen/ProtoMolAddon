/* -*- c++ -*- */
#ifndef SYSTEMFORCE_H
#define SYSTEMFORCE_H

#include <protomol/force/Force.h>

namespace ProtoMol {
  class ForceGroup;
  class GenericTopology;
  class Vector3DBlock;
  class ScalarStructure;
  class CompareForce;
  class TimeForce;
  //________________________________________ SystemForce

  class SystemForce : virtual public Force {
    // This class contains the definition of one system force (a force that 
    // works only on positions)

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    SystemForce() : Force() {}
    virtual ~SystemForce() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    virtual void evaluate(const GenericTopology *topo,
                          const Vector3DBlock *positions,
                          Vector3DBlock *forces,
                          ScalarStructure *energies); //=0;
    // Evaluate this force using the information given.

    virtual void evaluate(GenericTopology *topo,
                          const Vector3DBlock *positions,
                          Vector3DBlock *forces,
                          ScalarStructure *energies);

    virtual void parallelEvaluate(const GenericTopology *topo,
                                  const Vector3DBlock *positions,
                                  Vector3DBlock *forces,
                                  ScalarStructure *energies);

    //    virtual void computeBorn();
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

    virtual CompareForce *makeCompareForce(Force *actualForce,
                                           CompareForce *compareForce) const;
    virtual TimeForce *makeTimeForce(Force *actualForce) const;

    virtual void addToForceGroup(ForceGroup *forceGroup);
  };

  //________________________________________ INLINES

  inline void SystemForce::evaluate(const GenericTopology *topo,
                                    const Vector3DBlock *positions,
                                    Vector3DBlock *forces,
                                    ScalarStructure *energies) {}

  // Using dual cast to avoid code replication
  // for const and non-const
  // Ref. Scott Meyers, Effective C++ Third Edition
  // pp. 23-24 (C) 2007 Addison-Wesley
  inline void SystemForce::evaluate(GenericTopology *topo,
                                    const Vector3DBlock *positions,
                                    Vector3DBlock *forces,
                                    ScalarStructure *energies) {
    evaluate(static_cast<const GenericTopology *>(topo),
             positions, forces, energies);
  }

  inline void SystemForce::parallelEvaluate(const GenericTopology *topo,
                                            const Vector3DBlock *positions,
                                            Vector3DBlock *forces,
                                            ScalarStructure *energies) {
    evaluate(topo, positions, forces, energies);
  }
}
#endif /* SYSTEMFORCE_H */



