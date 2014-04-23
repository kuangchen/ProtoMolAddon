/*  -*- c++ -*-  */
#ifndef _LEAPFROG_SYMPATHETIC_COOLING_H
#define _LEAPFROG_SYMPATHETIC_COOLING_H

#include <protomol/addon/sympathetic_cooling/SympatheticCooling.h>
#include <protomol/addon/sympathetic_cooling/ThermalAtom.h>
#include <protomol/addon/sympathetic_cooling/IsotropicScattering.h>
#include <protomol/integrator/STSIntegrator.h>
#include <string>
#include <vector>

using namespace ProtoMolAddon::SympatheticCooling;
using ProtoMol::STSIntegrator;
using namespace std;


namespace ProtoMolAddon {
  //____ LeapfrogBufferGasIntegrator
  class LeapfrogSympatheticCoolingIntegrator : public STSIntegrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    LeapfrogSympatheticCoolingIntegrator();
    LeapfrogSympatheticCoolingIntegrator(Real timestep, const string& filename, ForceGroup *overloadedForces);

    ~LeapfrogSympatheticCoolingIntegrator();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LeapfrogSympatheticCoolingIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    void doKickdoDrift();
    void doHalfKickdoDrift();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return keyword;}
    virtual void getParameters(std::vector<Parameter> &parameters) const;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp *app);
    virtual void run(int numTimesteps);
    virtual void updateBeta(Real dt);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual STSIntegrator *doMake(const std::vector<Value> &values,
                                  ForceGroup *fg) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    SympatheticCooling::SympatheticCooling<ThermalAtom, IsotropicScattering> cooling;
    string fname;

  public:
    static const string keyword;
  };
}

#endif

