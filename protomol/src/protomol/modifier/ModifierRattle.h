/*  -*- c++ -*-  */
#ifndef MODIFIERRATTLE_H
#define MODIFIERRATTLE_H

#include <protomol/modifier/ModifierMetaRattle.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  class Integrator;

  //____ ModifierRattle
  class ModifierRattle : public ModifierMetaRattle {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierRattle();
    ModifierRattle(Real eps, int maxIter, bool all = true, int order = Constant::MAX_INT - 400);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return "Rattle";}
    virtual Modifier *doMake(const std::vector<Value> &values) const {
      return new ModifierRattle(values[0], values[1], values[2]);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual void doExecute(Integrator *i);
  };
}
#endif /* MODIFIERRATTLE_H */
