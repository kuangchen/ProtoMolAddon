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
#include <protomol/addon/LeapfrogBufferGasIntegrator2.h>
#include <protomol/addon/BufferGas.h>
#include <protomol/addon/Constants.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol::Constant;
using namespace ProtoMol;
using ProtoMolAddon::Constant::POSITION_CONV;

//____ LeapfrogBufferGasIntegrator2

const string LeapfrogBufferGasIntegrator2::keyword("LeapfrogBufferGas2");

LeapfrogBufferGasIntegrator2::LeapfrogBufferGasIntegrator2() :
  STSIntegrator() {}

LeapfrogBufferGasIntegrator2::LeapfrogBufferGasIntegrator2(Real timestep,
							   const string& filename, 
							   ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), 
  reader(filename),
  neutral_atom(reader),
  trap_radius(reader.GetValue<double>("trap.radius")),
  trap_z0(reader.GetValue<double>("trap.z0"))
{
  //std::cout << neutral_atom;
  //neutral_atom.SampleVelocityTest("test.log");
  //neutral_atom.ScatteringTest("test.log");
}

LeapfrogBufferGasIntegrator2::~LeapfrogBufferGasIntegrator2() {
}


void LeapfrogBufferGasIntegrator2::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
}

void LeapfrogBufferGasIntegrator2::doHalfKickdoDrift() {
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

void LeapfrogBufferGasIntegrator2::doKickdoDrift() {
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

void LeapfrogBufferGasIntegrator2::run(int numTimesteps) {
  double dt = getTimestep() / SI::TIME_FS;
  if (numTimesteps < 1)
    return;

  preStepModify();
  neutral_atom.CollideAll2(app, dt);
  //RemoveEnergeticIon(app);
  doHalfKickdoDrift();
  calculateForces();

  for (int i = 1; i < numTimesteps; i++) {
    neutral_atom.CollideAll2(app, dt);
    //RemoveEnergeticIon(app);
    doKickdoDrift();
    calculateForces();

  }
  doHalfKick();
  postStepModify();
}


STSIntegrator *LeapfrogBufferGasIntegrator2::doMake(const vector<Value> &values, ForceGroup *fg) const {
  return new LeapfrogBufferGasIntegrator2(values[0], values[1], fg);
}

//  --------------------------------------------------------------------  //
//  This function is necessary to compute the shadow Hamiltonian and it   //
//  is integrator specific.  This version is written to work with LF.     //
//  Update beta: beta -= dt * ( q * F + 2 U )                             //
//  --------------------------------------------------------------------  //

void LeapfrogBufferGasIntegrator2::updateBeta(Real dt) {
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

void LeapfrogBufferGasIntegrator2::getParameters(vector<Parameter> &parameters) const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back(Parameter("filename", Value(filename, ConstraintValueType::NotEmpty())));
}

void LeapfrogBufferGasIntegrator2::RemoveEnergeticIon(ProtoMolApp *app) {
  for (int j=0;j<app->positions.size();j++) {
    Vector3D r = app->positions[j] * POSITION_CONV;
    
    if ( (r[0]*r[0]+r[1]*r[1] > trap_radius * trap_radius || r[2]*r[2] > trap_z0 * trap_z0 ) 
	 && app->topology->atoms[j].scaledCharge > 0)
      app->topology->atoms[j].scaledCharge = 0;
    
  }
}
