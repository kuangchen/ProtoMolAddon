/*  -*- c++ -*-  */
#ifndef NPTVERLETINTEGRATOR_H
#define NPTVERLETINTEGRATOR_H

#include <protomol/integrator/STSIntegrator.h>

namespace ProtoMol {

  //_______________________________________________________________ NPTVerletIntegrator
  class GenericTopology;
  class ScalarStructure;
  class Vector3DBlock;
  class ForceGroup;
  template<class TIntegrator>
  class ModifierPreForceThermostat;
  template<class TIntegrator>
  class ModifierPostForceThermostat;
  template<class TIntegrator>
  class ModifierPreForceBarostat;
  template<class TIntegrator>
  class ModifierPostForceBarostat;
  class Modifer;

  class NPTVerletIntegrator: public STSIntegrator {
    template<class TIntegrator>
    friend class ModifierPreForceThermostat;
    template<class TIntegrator>
    friend class ModifierPostForceThermostat;
    template<class TIntegrator>
    friend class ModifierPreForceBarostat;
    template<class TIntegrator>
    friend class ModifierPostForceBarostat;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    NPTVerletIntegrator();  
    NPTVerletIntegrator(Real timestep,
			Real temperature,
			Real pressure,
			Real tauT,
			Real tauV,
			Real tauP,
			ForceGroup *overloadedForces);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class NPTVerletIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void do2ndHalfKick();
    void PreForceThermostat(); 
    void PostForceThermostat(); 
    void PreForceBarostat(); 
    void PostForceBarostat();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return 6;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Integrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void initialize(ProtoMolApp *);
    virtual void run(int numTimesteps);

    /// Create a Rattle modifier
    virtual Modifier* createRattleModifier(Real eps, int maxIter);
    /// Create a Shake modifier 
    virtual Modifier* createShakeModifier(Real eps, int maxIter);

  protected:
    virtual void addModifierAfterInitialize();
                  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class STSIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual STSIntegrator* doMake(const std::vector<Value>& values, ForceGroup* fg)const;
  protected:
    virtual void doDrift();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class StandardIntegrator
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    virtual void doHalfKick();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;
    Real getEpsilonVel() const {return myEpsilonVel;}
    Real getEtaVel()     const {return myEtaVel;}
    Real getNumAtoms()   const {return NumAtoms;}
    
  private:
    const Real myTargetTemp; //  Target temperature.  Units: (K)
    const Real myTargetPres; //  Target pressure.  Units: (bar)
    const Real myTauT;       //  thermostat oscillation period.  Units: (fs)
    const Real myTauV;       //  box volume thermostat oscillation period.  Units: (fs)
    const Real myTauP;       //  box volumr oscillation time constant.  Units: (fs)
    const Real kbT;          //  Target temperature multiplied by Boltzmann's constant.  Units: (kcal/mol)
    unsigned int NumAtoms;   //  Total # of atoms in the system.
    unsigned int myNumFree;  //  Total # of degrees of freedom = (3*Natoms - 3) - NumConstraints
    Real myVolume;           //  Current cubic volume.  Units (AA^3)
    Real myEpsilonVel;       //  Barostat strain rate velocity.  Units: (fs)^-1
    Real Qo;                 //  Particle thermostat mass.  Units: (kcal fs^2 / mol)
    Real Qv;                 //  Volume thermostat mass.  Units: (kcal fs^2 / mol)
    Real W;                  //  Barostat mass.  Units: (kcal fs^2 / mol)
    Real myEta;              //  Nose-Hoover particle thermostat variable.  Units: (dimensionless)
    Real myEtaV;             //  Nose-Hoover volume thermostat variable.  Units: (dimensionless)
    Real myEtaVel;           //  Velocity of the thermostat variable.  Units: (fs)^-1
    Real myEtaVolVel;        //  Velocity of the volume thermostat variable.  Units: (fs)^-1
  };

}
#endif

