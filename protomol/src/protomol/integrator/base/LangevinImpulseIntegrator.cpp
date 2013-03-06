#include <protomol/integrator/base/LangevinImpulseIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>

using namespace std; 
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ LangevinImpulseIntegrator

const string LangevinImpulseIntegrator::keyword("LangevinImpulse");

LangevinImpulseIntegrator::LangevinImpulseIntegrator() :
  STSIntegrator(), myLangevinTemperature(-1.0), myGamma(-1.0),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(-1) {}

LangevinImpulseIntegrator::
LangevinImpulseIntegrator(Real timestep, Real LangevinTemperature, Real gamma,
                          int seed, ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces),
  myLangevinTemperature(LangevinTemperature),
  myGamma(gamma / (1000 * Constant::INV_TIMEFACTOR)),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(seed) {}

void LangevinImpulseIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
}

void LangevinImpulseIntegrator::doDrift() {
  fluctuationLI();
  buildMolecularCenterOfMass(&app->positions, app->topology);
  buildMolecularMomentum(&app->velocities, app->topology);
}

// fluctuation using Dr. Skeel's LI scheme which involves a semi-update
// of velocities and a complete update of positions
void LangevinImpulseIntegrator::fluctuationLI() {
  const Real dt = getTimestep() * Constant::INV_TIMEFACTOR;   // in fs
  const Real tau1 = (1.0 - exp(-myGamma * dt)) / myGamma;
  const Real tau2 = (1.0 - exp(-2 * myGamma * dt)) / (2 * myGamma);
  const Real forceConstant =
    2 * Constant::BOLTZMANN * myLangevinTemperature * myGamma;
  const Real sqrtTau2 = sqrt(tau2);

  //  It is possible that roundoff and/or truncation can make this value < 0.
  //  It should be set to 0 if this happens.
  Real sqrtVal1 = (dt - tau1 * tau1 / tau2);

  if (sqrtVal1 < 0.) sqrtVal1 = 0;
  else sqrtVal1 = sqrt(sqrtVal1);

  for (unsigned int i = 0; i < app->positions.size(); i++) {
    Real mass = app->topology->atoms[i].scaledMass;
    Real sqrtFCoverM = sqrt(forceConstant / mass);
    Real langDriftVal = sqrtFCoverM / myGamma;
    Real langDriftZ1 = langDriftVal * (tau1 - tau2) / sqrtTau2;
    Real langDriftZ2 = langDriftVal * sqrtVal1;

    //  Generate gaussian random numbers for each spatial directions such that
    //  the average is 0 and the standard deviation is fParam.  The algorithm
    //  that is  used was taken from X-PLOR. The algorithm is a "sum of
    //  uniform deviates algorithm" which may be found in Abramowitz and
    //  Stegun, "Handbook of Mathematical Functions", pg 952.
    //force order of generation
    Real rand1 = randomGaussianNumber(mySeed);
    Real rand2 = randomGaussianNumber(mySeed);
    Real rand3 = randomGaussianNumber(mySeed);
      
    //into vector
    Vector3D gaussRandCoord1(rand3, rand2, rand1);
      
    rand1 = randomGaussianNumber(mySeed);
    rand2 = randomGaussianNumber(mySeed);
    rand3 = randomGaussianNumber(mySeed);
    Vector3D gaussRandCoord2(rand3, rand2, rand1);

    // update drift(fluctuation)
    app->positions[i] +=
      (gaussRandCoord1 * langDriftZ1 + gaussRandCoord2 * langDriftZ2 +
       app->velocities[i]) * tau1;

    // semi-update velocities
    app->velocities[i] = app->velocities[i] * exp(-myGamma * dt) +
      gaussRandCoord1 * sqrtFCoverM * sqrtTau2;
  }
}

void LangevinImpulseIntegrator::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("temperature", Value(myLangevinTemperature,
                                    ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("gamma", Value(myGamma * (1000 * Constant::INV_TIMEFACTOR),
                              ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("seed", Value(mySeed, ConstraintValueType::NotNegative()),
               1234));
}

STSIntegrator *LangevinImpulseIntegrator::doMake(const vector<Value> &values,
                                                 ForceGroup *fg) const {
  return new LangevinImpulseIntegrator(values[0], values[1], values[2],
                                       values[3], fg);
}
