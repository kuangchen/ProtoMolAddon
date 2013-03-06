/* -*- c++ -*- */
#ifndef INTEGRATOR_FACTORY_H
#define INTEGRATOR_FACTORY_H

#include <protomol/base/Factory.h>
#include <protomol/integrator/Integrator.h>

namespace ProtoMol {
  class ForceFactory;

  //________________________________________ IntegratorFactory
  class IntegratorFactory : public Factory<Integrator> {
    struct IntegratorInput {
      const Integrator *prototype;
      std::vector<Value> values;
      std::vector<std::string> forces;

      IntegratorInput(): prototype(0) {}
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    IntegratorFactory() {}
    virtual ~IntegratorFactory() {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From Factory<Integrator>
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    virtual void registerHelpText() const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class DetailsIntegratorFactory
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Integrator *make(const std::string &definition,
                     ForceFactory *forceFactory) const;
  };
}
#endif /* INTEGRATOR_FACTORY_H */
