/*  -*- c++ -*-  */
#ifndef GPU_H
#define GPU_H

#include <protomol/integrator/STSIntegrator.h>

#include <protomol/base/TimerStatistic.h>

namespace ProtoMol {
  class ScalarStructure;
  class ForceGroup;

  //____ GPU
  class GPU : public STSIntegrator {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    enum {MAXEQNN = 7};
    enum {NUMSW = 5};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    GPU();
    GPU(Real timestep, bool emu, int diag, ForceGroup *overloadedForces);
    ~GPU();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class GPU
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
    // New methods for class GPU
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    void doSteps(int numTimesteps);
    void doStepsE(int numTimesteps);
    void findForces(ForceGroup *overloadedForces);
    void gpuCalculateForces();
    Real minimalDist(Real *p1, Real *p2, Real *pbcCell, Real * pbcOrig, Real *mDiff);
    void outputDiagnostics(Real force, Real energy, Real value, Real deriv, Real distSquared, int i, int j);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;  
  private:
    //Diagnostic data
    Timer gpuTime;
    bool emulate;
    int diagnostics;
    //force field components    
    bool myBond;
    bool myAngle;
    bool myCoulomb;
    bool myCoulombDielec;
    bool myCoulombSCPISM;
    bool myCoulombBornRadii;
    bool myLennardJones;
    bool myDihedral;
    bool myImproper;
    Real cCutoff, cSwitchon, cSwitch, cOrder, cSwitchoff;   //switch data
    Real lCutoff, lSwitchon, lSwitch, lOrder, lSwitchoff;
    Real D, S, epsi;
    unsigned int sz; //size
    int swt;
    //Smooth switches
    static Real swcoef[][MAXEQNN], dswcoef[][MAXEQNN], d2swcoef[][MAXEQNN];
    int ordIdx;
    Real myIRange[MAXEQNN];
    //PBC
    Real *pbcCell, *hPbcCell, *pbcOrig, *mDiff, pbcMin2;
    bool myPeriodic;

  };
}

#endif

