#include <protomol/integrator/leapfrog/PLeapfrogIntegrator.h>
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

//____ PLeapfrogIntegrator

const string PLeapfrogIntegrator::keyword("PLeapfrog");

PLeapfrogIntegrator::PLeapfrogIntegrator() : STSIntegrator(), myTempForces(0) {}

PLeapfrogIntegrator::PLeapfrogIntegrator(Real timestep,
                                         ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces),
  myTempForces(new Vector3DBlock()) {}

PLeapfrogIntegrator::~PLeapfrogIntegrator() {
  if (myTempForces != 0) delete myTempForces;
}

void PLeapfrogIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
  myTempForces->resize(app->positions.size());
}

void PLeapfrogIntegrator::doKickdoDrift() {
  if (anyPreDriftOrNextModify() || anyPreStepModify() ||
      anyPostStepModify()) {
    doKick();

    if (anyPreStepModify() || anyPostStepModify()) {
      preDriftOrNextModify();
      doHalfDrift();
      postStepModify();
      preStepModify();
      doHalfDrift();
      postDriftOrNextModify();

    } else doDrift();
  } else {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] +=
        (*myForces)[i] * h / app->topology->atoms[i].scaledMass;
      app->positions[i] += app->velocities[i] * h;
    }

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
    postDriftOrNextModify();
  }
}

void PLeapfrogIntegrator::doKickdoHalfDrift() {
  if (anyPreDriftOrNextModify()) {
    doKick();
    preDriftOrNextModify();
    doHalfDrift();
  } else {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    for (unsigned int i = 0; i < count; ++i) {
      app->velocities[i] +=
        (*myForces)[i] * h / app->topology->atoms[i].scaledMass;
      app->positions[i] += app->velocities[i] * 0.5 * h;
    }

    buildMolecularCenterOfMass(&app->positions, app->topology);
    buildMolecularMomentum(&app->velocities, app->topology);
  }
}

void PLeapfrogIntegrator::doHalfDrift() {
  Real h = getTimestep() * Constant::INV_TIMEFACTOR;
  const unsigned int count = app->positions.size();

  for (unsigned int i = 0; i < count; ++i)
    app->positions[i] += app->velocities[i] * 0.5 * h;

  buildMolecularCenterOfMass(&app->positions, app->topology);
}

void PLeapfrogIntegrator::run(int numTimesteps) {
  if (numTimesteps < 1)
    return;

  preStepModify();
  doHalfDrift();
  postDriftOrNextModify();
  calculateForces();
  for (int i = 1; i < numTimesteps; ++i) {
    doKickdoDrift();
    calculateForces();
  }

  doKickdoHalfDrift();

  // Correction of energy ..
  Vector3DBlock *temp = myForces;
  myForces = myTempForces;
  Real t = app->topology->time;
  calculateForces();
  app->topology->time = t;
  myForces = temp;

  postStepModify();
}

STSIntegrator *PLeapfrogIntegrator::doMake(const vector<Value> &values,
                                           ForceGroup *fg) const {
  return new PLeapfrogIntegrator(values[0], fg);
}
