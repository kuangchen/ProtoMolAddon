/*  -*- c++ -*-  */
#ifndef MODIFIERREMOVELINEARMOMENTUM_H
#define MODIFIERREMOVELINEARMOMENTUM_H

#include <protomol/modifier/Modifier.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/ProtoMolApp.h>

namespace ProtoMol {
  //____ ModifierRemoveLinearMomentum
  class ModifierRemoveLinearMomentum : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierRemoveLinearMomentum(int freq = 0) :
      Modifier(Constant::MAX_INT - 100), myStep(0), myFreq(freq) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual void getParameters(std::vector<Parameter> &parameters) const {
      parameters.push_back
        (Parameter("-freq",
                   Value(myFreq, ConstraintValueType::NoConstraints())));
    }
    virtual std::string getIdNoAlias() const {return "RemoveLinearMomentum";}
    virtual Modifier *doMake(const std::vector<Value> &values) const {
      return new ModifierRemoveAngularMomentum(values[0]);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return false;}

  private:
    virtual void doExecute(Integrator *i) {
      if (myFreq == 0 || 0 == (myStep = (myStep % myFreq)))
        removeLinearMomentum(&app->velocities, app->topology);
      myStep++;
    }

  private:
    virtual void doInitialize() {myStep = 0;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    int myStep;
    int myFreq;
  };
}
#endif /* MODIFIER_H */
