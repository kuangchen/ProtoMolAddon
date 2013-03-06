/*  -*- c++ -*-  */
#ifndef MODIFIERINCREMENTTIMESTEP_H
#define MODIFIERINCREMENTTIMESTEP_H

#include <protomol/modifier/Modifier.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/ProtoMolApp.h>

namespace ProtoMol {
  class ModifierIncrementTimestep : public Modifier {
    static const std::string keyword;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierIncrementTimestep() : Modifier(Constant::MAX_INT) {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual void getParameters(std::vector<Parameter> &parameters) const {}
    virtual std::string getIdNoAlias() const {return "IncrementTimestep";}
    virtual Modifier *doMake(const std::vector<Value> &values) const {
      return new ModifierIncrementTimestep();
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual bool isInternal() const {return true;}

  private:
    virtual void doExecute(Integrator *i) {
      app->topology->time += i->getTimestep();
    }
  };
}
#endif /* MODIFIERINCREMENTTIMESTEP_H */
