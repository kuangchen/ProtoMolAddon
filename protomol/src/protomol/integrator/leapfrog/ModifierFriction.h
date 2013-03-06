/*  -*- c++ -*-  */
#ifndef MODIFIERFRICTION_H
#define MODIFIERFRICTION_H

#include <protomol/base/Exception.h>
#include <protomol/modifier/Modifier.h>
#include <protomol/integrator/leapfrog/NoseNVTLeapfrogIntegrator.h>

namespace ProtoMol {
  //____ ModifierFriction
  class ModifierFriction : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierFriction() : Modifier(Constant::MAX_INT - 400) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
    virtual std::string getIdNoAlias() const {
      return std::string("Friction");
    };

  private:
    virtual void doExecute(Integrator *i) {
      NoseNVTLeapfrogIntegrator *integrator = 
        dynamic_cast<NoseNVTLeapfrogIntegrator *>(i);
      if (!integrator) THROW("NoseNVTLeapfrogIntegrator required");
      integrator->friction();
    }
  };
}
#endif /* MODIFIER_H */
