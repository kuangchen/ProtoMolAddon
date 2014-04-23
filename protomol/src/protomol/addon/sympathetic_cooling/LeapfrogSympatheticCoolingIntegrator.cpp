#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/module/MainModule.h>
#include <algorithm>
#include <protomol/addon/sympathetic_cooling/LeapfrogSympatheticCoolingIntegrator.h>
#include <protomol/addon/sympathetic_cooling/SympatheticCooling.h>
#include <protomol/addon/Constants.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
using namespace ProtoMol::Constant;
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Constant;

//____ LeapfrogBufferGasIntegrator

const string LeapfrogSympatheticCoolingIntegrator::keyword("LeapfrogSympatheticCooling");

LeapfrogSympatheticCoolingIntegrator::LeapfrogSympatheticCoolingIntegrator() :
  STSIntegrator() {}

LeapfrogSympatheticCoolingIntegrator::LeapfrogSympatheticCoolingIntegrator(Real timestep,
							 const string& fname, 
							 ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), 
  cooling(fname, 0)
{
}

LeapfrogSympatheticCoolingIntegrator::~LeapfrogSympatheticCoolingIntegrator() {
}


void LeapfrogSympatheticCoolingIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();

  cooling.Initialize(app);
}

void LeapfrogSympatheticCoolingIntegrator::doHalfKickdoDrift() {
  if (anyPreDriftOrNextModify()) {
    doHalfKick();
    doDriftOrNextIntegrator();
  } else {
    Real h = getTimestep() * INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    //  Do a half kick on beta.
    updateBeta(0.5 * h);

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] += (*myForces)[i] * h * 0.5 /
                            app->topology->atoms[i].scaledMass;
      // app->positions[i] += app->velocities[i] * h;
    }
    
    app->positions += app->velocities*h;

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
    postDriftOrNextModify();
  }
}

void LeapfrogSympatheticCoolingIntegrator::doKickdoDrift() {
  if (anyPreDriftOrNextModify() || anyPreStepModify() ||
      anyPostStepModify()) {
    if (anyPreStepModify() || anyPostStepModify()) {
      doHalfKick();
      postStepModify();
      preStepModify();
      doHalfKick();
    } else
      doKick();
    doDriftOrNextIntegrator();
  } else {
    Real h = getTimestep() * INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    updateBeta(h);

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] +=
        (*myForces)[i] * h / app->topology->atoms[i].scaledMass;
      //  app->positions[i] += app->velocities[i] * h;
    }

    app->positions += app->velocities*h;

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
    postDriftOrNextModify();
  }
}

void LeapfrogSympatheticCoolingIntegrator::run(int numTimesteps) {
  double dt = getTimestep() / SI::TIME_FS;
  if (numTimesteps < 1)
    return;

  preStepModify();
  cooling.Run(dt);
  doHalfKickdoDrift();
  calculateForces();

  for (int i = 1; i < numTimesteps; i++) {
    cooling.Run(dt);
    doKickdoDrift();
    calculateForces();
  }

  doHalfKick();
  postStepModify();
}

STSIntegrator *LeapfrogSympatheticCoolingIntegrator::doMake(const vector<Value> &values, ForceGroup *fg) const {
  return new LeapfrogSympatheticCoolingIntegrator(values[0], values[1], fg);
}

//  --------------------------------------------------------------------  //
//  This function is necessary to compute the shadow Hamiltonian and it   //
//  is integrator specific.  This version is written to work with LF.     //
//  Update beta: beta -= dt * ( q * F + 2 U )                             //
//  --------------------------------------------------------------------  //

void LeapfrogSympatheticCoolingIntegrator::updateBeta(Real dt) {
  //  ----------------------------------------------------------------  //
  //  The shadow calculation is done in a postStep modifier.  If there  //
  //  aren't any, then obviously we don't need to do this calculation.  //
  //  It's possible that a different poststep modifier could make this  //
  //  execute, but no harm would be done ... only some extra cycles.    //
  //  ----------------------------------------------------------------  //

  if (!(anyPostStepModify() || top()->anyPostStepModify()))
    return;

  Real posDotF = 0.;

  for (unsigned int i = 0; i < app->positions.size(); i++)
    posDotF += app->positions[i].dot((*myForces)[i]);

  myBeta -= dt * (posDotF + 2. * myPotEnergy);
}

void LeapfrogSympatheticCoolingIntegrator::getParameters(vector<Parameter> &parameters) const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back(Parameter("filename", Value(fname, ConstraintValueType::NotEmpty())));
}


