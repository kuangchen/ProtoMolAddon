/*  -*- c++ -*-  */
#ifndef MODIFIERFORCEPROJECTION_H
#define MODIFIERFORCEPROJECTION_H

#include <protomol/modifier/Modifier.h>
#include <protomol/integrator/normal/NormalModeUtilities.h>

namespace ProtoMol {
  //____ ModifierForceProjection
  class ModifierForceProjection : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierForceProjection(NormalModeUtilities *i) :
      Modifier(Constant::MAX_INT - 400), myTheIntegrator(i) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
    virtual std::string getIdNoAlias() const {
      return std::string("ForceProjection");
    };

  private:
    virtual void doExecute(Integrator *i) {
      myTheIntegrator->forceProjection();
    }
    

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    NormalModeUtilities *myTheIntegrator;
  };
}
#endif /* MODIFIER_H */
