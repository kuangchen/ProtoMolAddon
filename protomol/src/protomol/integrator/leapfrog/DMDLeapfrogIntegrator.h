/*  -*- c++ -*-  */
#ifndef DMDLEAPFROGINTEGRATOR_H
#define DMDLEAPFROGINTEGRATOR_H

#include <protomol/integrator/STSIntegrator.h>

namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;
  class Vector3DBlock;

  /**
      For Dissipative MD, and it is self-consistent
   */
  //____DMDLeapfrogIntegrator
  class DMDLeapfrogIntegrator : public STSIntegrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    DMDLeapfrogIntegrator();
    /// Specify the the forces to evaluate.
    /// Also sigma, the noise level control parameter, a constant
    /// and numIter, the count of iterations to be made in the
    /// "self-consistency" loop. The desired temperature in Kelvin
    /// is used to get the drag coefficient, gamma.
    /// sigma^2 = 2 * gamma * BOLZMAN * Temperature.
    /// the gamma is passed in assuming unit of (ps^-1).
    /// Internally, it is used in the unit of (fs^-1), a factor of 0.001
    DMDLeapfrogIntegrator(Real timestep, int numIter, Real gamma,
                          Real initialTemperature, int seed,
                          ForceGroup *overloadedForces);
    ~DMDLeapfrogIntegrator();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class DMDLeapfrogIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void calculateDissipativeForces();
    void calculateDissipativeAndRandomForces();
    void doHalfKickVhat();
    void doHalfKickIterate();

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

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual void doHalfKick();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual STSIntegrator *doMake(const std::vector<Value> &values,
                                  ForceGroup *fg) const;

  protected:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
  private:
    Vector3DBlock *myDissipativeForces;
    Vector3DBlock *myRandomForces;
    Vector3DBlock *myVhat;
    Real myDissipativeCutoff;
    Real myGamma;
    Real myTemperature;
    int myNumIter;
    Real mySigma;
    int mySeed;
  };
  //____ INLINES
}
#endif

