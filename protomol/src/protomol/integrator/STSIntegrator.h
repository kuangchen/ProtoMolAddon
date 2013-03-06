/*  -*- c++ -*-  */
#ifndef STSINTEGRATOR_H
#define STSINTEGRATOR_H

#include <protomol/integrator/StandardIntegrator.h>
#include <protomol/topology/GenericTopology.h>

namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;

  //________________________________________ STSIntegrator

  class STSIntegrator : public StandardIntegrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    STSIntegrator();
    STSIntegrator(Real timestep, ForceGroup *overloadedForces);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    STSIntegrator *make(const std::vector<Value> &values, ForceGroup *fg) const;
    virtual Real setTimestep(Real);
  protected:
    virtual void doDrift();
  private:
    virtual STSIntegrator *doMake(const std::vector<Value> &values,
				  ForceGroup *fg) const = 0;


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void getParameters(std::vector<Parameter> &parameter) const;
    virtual std::string getIdNoAlias() const {return "STSIntegrator";}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual Integrator *next() {return NULL;}
    virtual const Integrator *next() const {return NULL;}
    virtual Real getTimestep() const;
  protected:
    virtual void addModifierAfterInitialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual void doDriftOrNextIntegrator();
    virtual void calculateForces();
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myTimestep;
  };
  //________________________________________ INLINES

  inline Real STSIntegrator::getTimestep() const {
    return isForward() ? myTimestep : -myTimestep;
  }

  inline Real STSIntegrator::setTimestep(Real newTimestep) {
    Real oldTimestep = myTimestep;

    myTimestep = newTimestep;

    return oldTimestep;
  }
}
#endif
