/*  -*- c++ -*-  */
#ifndef LEAPFROGBUFFERGASINTEGRATOR2_H
#define LEAPFROGBUFFERGASINTEGRATOR2_H

#include <protomol/integrator/STSIntegrator.h>
//#include <protomol/addon/BufferGasManager.h>
#include <protomol/addon/LuaConfigReader.h>
#include <protomol/addon/NeutralAtom.h>

#include <string>
#include <vector>
                
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Lua;
using ProtoMol::STSIntegrator;
using namespace std;


namespace ProtoMolAddon {
  //____ LeapfrogBufferGasIntegrator2
  class LeapfrogBufferGasIntegrator2 : public STSIntegrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    LeapfrogBufferGasIntegrator2();
    LeapfrogBufferGasIntegrator2(Real timestep, const string& filename, ForceGroup *overloadedForces);

    ~LeapfrogBufferGasIntegrator2();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class LeapfrogBufferGasIntegrator2
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
    string filename; 
    LuaConfigReader reader;
    NeutralAtom neutral_atom;
    double trap_radius;
    double trap_z0;
  public:
    static const string keyword;
  };
}

#endif

