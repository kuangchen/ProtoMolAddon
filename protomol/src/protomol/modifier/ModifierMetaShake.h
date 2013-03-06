/*  -*- c++ -*-  */
#ifndef MODIFIERMETASHAKE_H
#define MODIFIERMETASHAKE_H

#include <protomol/modifier/ModifierMetaRattleShake.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  class Integrator;
  //____ ModifierMetaShake
  class ModifierMetaShake : public ModifierMetaRattleShake {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierMetaShake(Real eps, int maxIter, bool all, int order);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class ModifierMetaShakeShake
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual Real calcError() const;
  };
}
#endif /* MODIFIERMETASHAKE_H */
