/*  -*- c++ -*-  */
#ifndef MODIFIERSHAKE_H
#define MODIFIERSHAKE_H

#include <protomol/modifier/ModifierMetaShake.h>
#include <protomol/type/Vector3DBlock.h>

namespace ProtoMol {
  class Integrator;

  //____ ModifierShake
  class ModifierShake : public ModifierMetaShake {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierShake();
    ModifierShake(Real eps, int maxIter, bool all = true, int order = Constant::MAX_INT - 400);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual std::string getIdNoAlias() const {return "Shake";}
    virtual Modifier *doMake(const std::vector<Value> &values) const {
      return new ModifierShake(values[0], values[1], values[2]);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual void doExecute(Integrator *i);
  };
}
#endif /* MODIFIERSHAKE_H */
