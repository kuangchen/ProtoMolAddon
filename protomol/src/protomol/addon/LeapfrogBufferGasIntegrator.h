/*  -*- c++ -*-  */
#ifndef LEAPFROGBUFFERGASINTEGRATOR_H
#define LEAPFROGBUFFERGASINTEGRATOR_H

#include <protomol/integrator/STSIntegrator.h>
#include <protomol/addon/BufferGasManager.h>
#include <protomol/addon/LuaConfigReader.h>
#include <string>
#include <vector>
                
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Lua;
using ProtoMol::STSIntegrator;
using namespace std;


namespace ProtoMolAddon {
  //____ LeapfrogBufferGasIntegrator
  class LeapfrogBufferGasIntegrator : public STSIntegrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    LeapfrogBufferGasIntegrator();
    LeapfrogBufferGasIntegrator(Real timestep, const string& filename, ForceGroup *overloadedForces);

    ~LeapfrogBufferGasIntegrator();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LeapfrogBufferGasIntegrator
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
    void RemoveEnergeticIon(ProtoMolApp *app);
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
    BufferGasManager bg_manager;
    string filename; 
    LuaConfigReader reader;
    double trap_radius;
  public:
    static const string keyword;
    void runGeneral(Real t1, Real t2);
  };
}

#endif

