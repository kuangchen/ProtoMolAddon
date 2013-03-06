/*  -*- c++ -*-  */
#ifndef MODIFIERMETARATTLE_H
#define MODIFIERMETARATTLE_H

#include <protomol/modifier/ModifierMetaRattleShake.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  class Integrator;
  //____ ModifierMetaRattle
  class ModifierMetaRattle : public ModifierMetaRattleShake {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierMetaRattle(Real eps, int maxIter, bool all, int order);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class ModifierMetaRattleShake
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual Real calcError() const;
  };
}
#endif /* MODIFIERMETARATTLE_H */
