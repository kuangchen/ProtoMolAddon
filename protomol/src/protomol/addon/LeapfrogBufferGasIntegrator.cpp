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
#include <vector>
#include <protomol/addon/LeapfrogBufferGasIntegrator.h>
#include <protomol/addon/BufferGas.h>
#include <protomol/addon/Constants.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
using namespace ProtoMol::Constant;
using namespace ProtoMolAddon;
using namespace ProtoMolAddon::Constant;

//____ LeapfrogBufferGasIntegrator

const string LeapfrogBufferGasIntegrator::keyword("LeapfrogBufferGas");

LeapfrogBufferGasIntegrator::LeapfrogBufferGasIntegrator() :
  STSIntegrator() {}

LeapfrogBufferGasIntegrator::LeapfrogBufferGasIntegrator(Real timestep,
							 const string& filename, 
							 ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), 
  bg_manager(),
  reader(filename),
  trap_radius(reader.GetValue<double>("trap.radius"))
{
}

LeapfrogBufferGasIntegrator::~LeapfrogBufferGasIntegrator() {
}


void LeapfrogBufferGasIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();

  bg_manager.InitializeBufferGas(reader);
  bg_manager.InitializeCollisionSchedule(reader, app);
}

void LeapfrogBufferGasIntegrator::doHalfKickdoDrift() {
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


// Run simulations n steps and eta fraction step
void LeapfrogBufferGasIntegrator::runGeneral(double t1, double t2) {
  //std::cout << "Running from " << t1 << " to " << t2 << "\n";

  int n;
  Real eta;

  Real stepsize = getTimestep() / SI::TIME_FS;
  n = static_cast<int> ((t2-t1)/stepsize + 0.5);
  eta = (t2-t1)/stepsize - n;

  if (n>0) {
    doHalfKickdoDrift();
    RemoveEnergeticIon(app);
    calculateForces();
    for (int i=1; i<n; i++) {
      doKickdoDrift();
      RemoveEnergeticIon(app);
      calculateForces();
    }
    doHalfKick();
  }

  if (fabs(eta) > EPSILON) {
    double oldTimestep = getTimestep();
    setTimestep(oldTimestep * eta);
  
    doHalfKickdoDrift ();
    RemoveEnergeticIon(app);
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

  double stepSize = getTimestep() / SI::TIME_FS;
  double current = app->currentStep * stepSize;
  double end = (app->currentStep + numTimesteps) * stepSize;

  double next_collision_time = 0;

  while (!bg_manager.IsCollisionFinished() && 
	 (next_collision_time = bg_manager.GetNextCollisionTime()) < end) {

    runGeneral(current, next_collision_time);
    bg_manager.Collide(app);
    current = next_collision_time;
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

void LeapfrogBufferGasIntegrator::getParameters(vector<Parameter> &parameters) const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back(Parameter("filename", Value(filename, ConstraintValueType::NotEmpty())));
}

void LeapfrogBufferGasIntegrator::RemoveEnergeticIon(ProtoMolApp *app) {
  for (int j=0;j<app->positions.size();j++) {
    Vector3D r = app->positions[j] * POSITION_CONV;
    
    if (r[0]*r[0]+r[1]*r[1] > trap_radius * trap_radius && app->topology->atoms[j].scaledCharge > 0) 
      app->topology->atoms[j].scaledCharge = 0;
    
  }
}

