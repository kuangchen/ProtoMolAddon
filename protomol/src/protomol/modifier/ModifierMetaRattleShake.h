/*  -*- c++ -*-  */
#ifndef MODIFIERMETARATTLESHAKE_H
#define MODIFIERMETARATTLESHAKE_H

#include <protomol/modifier/Modifier.h>
#include <protomol/topology/Bond.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  //____ ModifierMetaRattleShake

  /**
      Base class of Rattle and Shake algorithms. It implements the construction
      of the constraints, provides the maximal common interface and data 
      members. The meta classes of Rattle and Shake implement the error
      estimate, where the actual constraints are implemented in the concrete
      class (for NVE) or used by a concrete templated class (NPT, since
      getEpsilonVel() and getEtaVel() must be accessible). myListOfConstraints
      is pointer to the actual list.
   */
  class ModifierMetaRattleShake : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierMetaRattleShake(Real eps, int maxIter, bool all, int order);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ModifierMetaRattleShake
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual Real calcError() const = 0;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return false;}
    virtual void getParameters(std::vector<Parameter> &parameters) const;

  private:
    virtual void doInitialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    Real myEpsilon;
    int myMaxIter;
    bool myAll;

    Vector3DBlock myLastPositions;
    const std::vector<Bond::Constraint> *myListOfConstraints;
  };
}
#endif /* MODIFIERMETARATTLESHAKE_H */
