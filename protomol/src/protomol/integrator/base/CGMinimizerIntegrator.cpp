#include <protomol/integrator/base/CGMinimizerIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ CGMinimizerIntegrator

const string CGMinimizerIntegrator::keyword("CGMinimizer");

CGMinimizerIntegrator::CGMinimizerIntegrator() :
  STSIntegrator() {
  gk = 0; pk = 0; oldPos = 0;
}

CGMinimizerIntegrator::CGMinimizerIntegrator(Real timestep, Real alpha, 
                                             Real beta2, Real restart,
                                             ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), myAlpha(alpha), myBeta(beta2),
  myRestart(restart) {
  gk = 0; pk = 0; oldPos = 0;
}

CGMinimizerIntegrator::~CGMinimizerIntegrator() {
  if (gk != 0) delete gk;
  if (pk != 0) delete pk;
  if (oldPos != 0) delete oldPos;
}

void CGMinimizerIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
  //CG initialization
  gk = new Vector3DBlock;
  pk = new Vector3DBlock;
  oldPos = new Vector3DBlock;
  calculateForces();
  (*gk).resize(app->positions.size());
  (*gk).intoSubtract(*myForces);
  (*pk).resize(app->positions.size());
  (*pk).intoSubtract(*gk);
  beta = 0;
  lastLambd = 1.0;
  mcount = 0;
  //
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Non-linear Conjugate Gradient Minimizer.
// Uses the Fletcher-Reeves Algarithm for \beta
// The line-search section works by using the Newton-Raphson
// method to find the root of g_{k+1}^Tp_k which is the slope
// of the potential energy w.r.t. \lambda, U(\lambda).
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void CGMinimizerIntegrator::run(int numTimesteps) {
  Real b1, b2, U1, dU, gp, gp1p, gpr, oldGp1p, tempf, oldLambda;
  int i, k;

  if (numTimesteps < 1)
    return;

  for (int j = 0; j < numTimesteps; j++) {
    mcount++;
    //calculate part of beta before gk changes
    b1 = b2 = 0.0;
    for (i = 0; i < (int)app->positions.size(); i++)
      b2 += (*gk)[i].dot((*gk)[i]);

    //do line search for lambda
    U1 = app->energies.potentialEnergy();
    gp = 0.0;
    for (i = 0; i < (int)app->positions.size(); i++)
      gp += (*gk)[i].dot((*pk)[i]);

    gp1p = gp;
    oldLambda = 0.0;
    //calc initial lambda
    lambda = lastLambd;   //0.0001;//-1.0/(gp*myAlpha);
    //save positions
    *oldPos = app->positions;
    //
    gpr = dU = 2.0;   //preset to pass while test
    // meet the |g_{k+1}^Tp_k|=<|g_k^Tp_k|
    k = 0;
    do {
      if (++k > 8) break;
      app->positions = *oldPos;
      app->positions.intoWeightedAdd(lambda, *pk);
      calculateForces();
      dU = app->energies.potentialEnergy() - U1;
      (*gk).intoWeighted(-1.0, *myForces);
      oldGp1p = gp1p;
      gp1p = 0.0;
      for (i = 0; i < (int)app->positions.size(); i++)
        gp1p += (*gk)[i].dot((*pk)[i]);

      gpr = fabs(gp1p / gp);
      if (gpr > myBeta) {
        tempf = lambda;
        lambda = lambda - ((lambda - oldLambda) / (gp1p - oldGp1p)) * gp1p;
        oldLambda = tempf;
      }
    } while (gpr > myBeta);

    //
    if (k > 8 || dU > (myAlpha * lambda * gp) ||
        (myRestart > 0 && mcount > (int)myRestart)) {                           
      // re-start method if stuck
      j = numTimesteps;
      app->positions = *oldPos;
      calculateForces();
      (*gk).intoWeighted(-1.0, *myForces);
      //update pk
      (*pk).intoWeighted(-1.0, *gk);
      mcount = 0;
      lastLambd = 1.0;
      oldLambda = 0.0;
      oldGp1p = gp;
    } else {   // or do normal p_k update
      for (i = 0; i < (int)app->positions.size(); i++)
        b1 += (*gk)[i].dot((*gk)[i]);

      if (b2) beta = b1 / b2;
      else beta = 0.0;
      if (mcount == 1) lastLambd = lambda;
      if (lambda > lastLambd) lastLambd = lambda;
      //update pk
      (*pk).intoWeighted(beta, *pk);
      (*pk).intoSubtract(*gk);
    }
  }
}

void CGMinimizerIntegrator::getParameters(vector<Parameter> &parameters) const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("alpha", Value(myAlpha, ConstraintValueType::NotNegative()),
               0.001,
               Text("steplength 'sufficient decrease' parameter")));
  parameters.push_back
    (Parameter("beta", Value(myBeta, ConstraintValueType::NotNegative()), 0.05,
               Text("steplength 'sufficient directional derivative reduction' "
                    "parameter")));
  parameters.push_back
    (Parameter("restart", Value(myBeta, ConstraintValueType::NotNegative()), 0,
               Text("steplength 'restart interval")));
}

STSIntegrator *CGMinimizerIntegrator::doMake(const vector<Value> &values,
                                             ForceGroup *fg) const {
  return new CGMinimizerIntegrator(values[0], values[1], values[2], values[3],
                                   fg);
}
