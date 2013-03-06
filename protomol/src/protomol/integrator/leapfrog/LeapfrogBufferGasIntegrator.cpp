#include <protomol/integrator/leapfrog/LeapfrogBufferGasIntegrator.h>
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
#include <cmath>
#include <string>
#include <protomol/integrator/leapfrog/BufferGas.h>

using std::string;
using std::sort;

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ LeapfrogBufferGasIntegrator

const string LeapfrogBufferGasIntegrator::keyword("LeapfrogBufferGas");

LeapfrogBufferGasIntegrator::LeapfrogBufferGasIntegrator() :
  STSIntegrator() {}

LeapfrogBufferGasIntegrator::LeapfrogBufferGasIntegrator(Real timestep,
							 const string& filename, 
							 ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), 
  L(filename),
  gas(L),
  nextEvent(0)
{
}

LeapfrogBufferGasIntegrator::~LeapfrogBufferGasIntegrator() {
}


bool compareFunction(pair<int, Real> p1, pair<int, Real> p2) {
  return p1.second < p2.second;
}

void LeapfrogBufferGasIntegrator::collide() {
  Real p_conv = 1e-10;
  Real v_conv = 1e-10 * ProtoMol::Constant::SI::TIME_FS / ProtoMol::Constant::TIMEFACTOR;
  int atomIndex = schedule[nextEvent].atom;
  Real m = app->topology->atoms[atomIndex].scaledMass * Constant::SI::AMU;
  Vector3D p = app->positions[atomIndex] * p_conv;
  Vector3D v = app->velocities[atomIndex] * v_conv; 

  gas.collide(m, p, v);
  app->velocities[atomIndex] = v / v_conv;
}


void LeapfrogBufferGasIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();

  Real step = getTimestep() / Constant::SI::TIME_FS;
  int atomCount = app->topology->atoms.size();
  Real start = (int)app->config[InputFirststep::keyword] * step;
  Real end = app->lastStep * step;

  schedule = gas.createCollisionSchedule(start, end, atomCount);
}

void LeapfrogBufferGasIntegrator::doHalfKickdoDrift() {
  if (anyPreDriftOrNextModify()) {
    doHalfKick();
    doDriftOrNextIntegrator();
  } else {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    const unsigned int count = app->positions.size();

    //  Do a half kick on beta.
    updateBeta(0.5 * h);

    #pragma omp parallel for
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

void LeapfrogBufferGasIntegrator::doKickdoDrift() {
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
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
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


// Run simulations n steps and eta fraction step
void LeapfrogBufferGasIntegrator::runGeneral(double t1, double t2) {
  int n;
  Real eta;

  Real stepsize = getTimestep() / ProtoMol::Constant::SI::TIME_FS;
  n = static_cast<int> ((t2-t1)/stepsize + 0.5);
  eta = (t2-t1)/stepsize - n;

  if (n>0) {
    doHalfKickdoDrift();
    calculateForces();
    for (int i=1; i<n; i++) {
      doKickdoDrift();
      calculateForces();
    }
    doHalfKick();
  }

  if (eta > Constant::EPSILON) {
    double oldTimestep = getTimestep();
    setTimestep(oldTimestep * eta);
  
    doHalfKickdoDrift ();
    calculateForces ();
    doHalfKick();
    calculateForces ();

    setTimestep (oldTimestep);
  }
  return;
}

void LeapfrogBufferGasIntegrator::run(int numTimesteps) {
  if (numTimesteps < 1)
    return;

  preStepModify();
  Real stepSize = getTimestep() / Constant::SI::TIME_FS;

  Real current = app->currentStep * stepSize;
  Real end = (app->currentStep + numTimesteps) * stepSize;
  Real next = 0;

  while (nextEvent < schedule.size() && (next = schedule[nextEvent].time) < end) {
    runGeneral(current, next);
    current = next;
    collide();
    nextEvent++;
  }

  runGeneral(current, end);
  postStepModify();
}

STSIntegrator *LeapfrogBufferGasIntegrator::doMake(const vector<Value> &values, ForceGroup *fg) const {
  return new LeapfrogBufferGasIntegrator(values[0], values[1], fg);
}

//  --------------------------------------------------------------------  //
//  This function is necessary to compute the shadow Hamiltonian and it   //
//  is integrator specific.  This version is written to work with LF.     //
//  Update beta: beta -= dt * ( q * F + 2 U )                             //
//  --------------------------------------------------------------------  //

void LeapfrogBufferGasIntegrator::updateBeta(Real dt) {
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

void LeapfrogBufferGasIntegrator::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back(Parameter("filename", Value(filename, ConstraintValueType::NotEmpty())));
}
