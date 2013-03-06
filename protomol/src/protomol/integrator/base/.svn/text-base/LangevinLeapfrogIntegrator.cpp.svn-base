#include <protomol/integrator/base/LangevinLeapfrogIntegrator.h>
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
//____ LangevinLeapfrogIntegrator

const string LangevinLeapfrogIntegrator::keyword("LangevinLeapfrog");

LangevinLeapfrogIntegrator::LangevinLeapfrogIntegrator() :
  STSIntegrator(), myLangevinTemperature(-1.0), myGamma(-1.0),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(-1) {}

LangevinLeapfrogIntegrator::
LangevinLeapfrogIntegrator(Real timestep, Real LangevinTemperature, Real gamma,
                          int seed, ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces),
  myLangevinTemperature(LangevinTemperature),
  myGamma(gamma / (1000 * Constant::INV_TIMEFACTOR)),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(seed) {}

void LangevinLeapfrogIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
}

void LangevinLeapfrogIntegrator::doDrift() {
  const Real h = getTimestep() * Constant::INV_TIMEFACTOR;
  app->positions.intoWeightedAdd(h, app->velocities);
  buildMolecularCenterOfMass(&app->positions, app->topology);
  buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinLeapfrogIntegrator::doHalfKick() {
    const unsigned int count = app->positions.size();
    const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
    const Real fdt = ( 1.0 - exp( -0.5 * myGamma * dt ) ) / myGamma;
    const Real vdt = exp(-0.5*myGamma*dt);
    const Real ndt = sqrt( ( 1.0 - exp( -myGamma * dt ) ) / (2.0 * myGamma) );
    const Real forceConstant = 2 * Constant::BOLTZMANN * myLangevinTemperature *
      myGamma;

    for (unsigned int i = 0; i < count; i++ ) {
        //  Generate gaussian random numbers for each spatial direction
        //force order of generation
        Real rand1 = randomGaussianNumber(mySeed);
        Real rand2 = randomGaussianNumber(mySeed);
        Real rand3 = randomGaussianNumber(mySeed);
        
        //into vector
        Vector3D gaussRandCoord1(rand3, rand2, rand1);
        
        Real mass = app->topology->atoms[i].scaledMass;
        Real sqrtFCoverM = sqrt(forceConstant / mass);
        // semi-update velocities
        app->velocities[i] = app->velocities[i]*vdt
                                +(*myForces)[i] * fdt / mass
                                    +gaussRandCoord1*sqrtFCoverM*ndt;
    }
    buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinLeapfrogIntegrator::getParameters(vector<Parameter> &parameters)
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

STSIntegrator *LangevinLeapfrogIntegrator::doMake(const vector<Value> &values,
                                                 ForceGroup *fg) const {
  return new LangevinLeapfrogIntegrator(values[0], values[1], values[2],
                                       values[3], fg);
}
