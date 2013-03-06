/*  -*- c++ -*-  */
#ifndef NONSTANDARDINTEGRATOR_H
#define NONSTANDARDINTEGRATOR_H

#include <protomol/integrator/Integrator.h>

namespace ProtoMol {
  //________________________________________ NonStandardIntegrator
  class ForceGroup;

  class NonStandardIntegrator : public Integrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NonStandardIntegrator() {};
    NonStandardIntegrator(ForceGroup *forceGroup) : Integrator(forceGroup) {};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NonStandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual NonStandardIntegrator *make(std::vector<Value> values,
                                        ForceGroup *fg) const = 0;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
  };
  //________________________________________ INLINES
}
#endif
