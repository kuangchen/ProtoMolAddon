/*  -*- c++ -*-  */
#ifndef RMTINTEGRATOR_H
#define RMTINTEGRATOR_H

#include <protomol/integrator/STSIntegrator.h>
#include <iostream>
#include <stdio.h>

namespace ProtoMol {


  //_________________________________________________________________ RMTIntegrator

  class ScalarStructure;
  class ForceGroup;
  class Vector3DBlock;

  class RMTIntegrator : public STSIntegrator {

  private:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {MAX_THERMOSTAT = 5};


    //friend class ModifierFriction;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    RMTIntegrator();
    RMTIntegrator(Real timestep,
            Real temp,
            Real q1, Real q2, Real q3, Real q4, Real q5, 
            Real c2, Real c3, Real c4, Real c5, 
            int numst, bool tdof,
            ForceGroup *overloadedForces);
    ~RMTIntegrator();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class RMTIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
      virtual void streamRead( std::istream& inStream );
      virtual void streamWrite( std::ostream& outStream ) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp *app);
    virtual void run(int numTimesteps);

  protected:
    void halfUpdtH2(int typ);
    Real prodSs(int start, int end);
    void halfUpdtH3j(int dir);
    void halfUpdtH31();
    void UpdtH1();
    Real totalEnergy(int typ);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual STSIntegrator *doMake(const std::vector<Value> &values,
                                  ForceGroup *fg) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    const Real   myTemp;
    const Real   myQ1, myQ2, myQ3, myQ4, myQ5;
    const Real   myC2, myC3, myC4, myC5;
    int        myNumStats;
    Real       myQ[MAX_THERMOSTAT];
    Real       myC[MAX_THERMOSTAT];
    Real       myS[MAX_THERMOSTAT];
    Real       myPs[MAX_THERMOSTAT];
    Real		   myAvTKE[MAX_THERMOSTAT];
    Real		   myAvS[MAX_THERMOSTAT];
    Real		   myh0;
    Real		   myOldProdS;
    Real		   mySn, myNf, mykT;
    std::string	   tstatfile;
    int			   totStep;
    Real	     avKE, avKEsq;
    bool		   incTdof;

  };
  //______________________________________________________________________ INLINES
}

#endif
