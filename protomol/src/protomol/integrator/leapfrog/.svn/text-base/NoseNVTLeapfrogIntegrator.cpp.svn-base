#include <protomol/integrator/leapfrog/NoseNVTLeapfrogIntegrator.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/integrator/leapfrog/ModifierFriction.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ NoseNVTLeapfrogIntegrator

const string NoseNVTLeapfrogIntegrator::keyword("NoseNVTLeapfrog");

NoseNVTLeapfrogIntegrator::NoseNVTLeapfrogIntegrator() :
  STSIntegrator(), myTemperature(0.0), myThermalInertia(0.0),
  myBathPosition(0.0) {}

NoseNVTLeapfrogIntegrator::NoseNVTLeapfrogIntegrator(
  Real timestep, Real temperature, Real thermalInertia, Real bathPosition,
  ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), myTemperature(temperature),
  myThermalInertia(thermalInertia), myBathPosition(bathPosition) {}

void NoseNVTLeapfrogIntegrator::friction() {
  const unsigned int numberOfAtoms = app->topology->atoms.size();
  const Real h = getTimestep() * Constant::INV_TIMEFACTOR;
  Real kineticEnergy = 0.0;
  for (unsigned int i = 0; i < numberOfAtoms; ++i) {
    Real m = app->topology->atoms[i].scaledMass;
    //+(*myForces)[i]*h/m
    kineticEnergy += (app->velocities[i]).normSquared() * m;
  }

  myBathPosition +=
    (0.5 * kineticEnergy -
     myTargetKE) * myThermalInertia / (h * numberOfAtoms);

  for (unsigned int i = 0; i < numberOfAtoms; i++)
    (*myForces)[i] -= app->velocities[i] * myBathPosition *
                      app->topology->atoms[i].scaledMass;
}

void NoseNVTLeapfrogIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);

  myTargetKE =
    (3.0 / 2.0 * app->positions.size() * Constant::BOLTZMANN * myTemperature);
  mySumMass = 0.0;
  for (unsigned int i = 0; i < app->positions.size(); i++)
    mySumMass += app->topology->atoms[i].scaledMass;

  initializeForces();
}

void NoseNVTLeapfrogIntegrator::addModifierAfterInitialize() {
  adoptPostForceModifier(new ModifierFriction());
  STSIntegrator::addModifierAfterInitialize();
}

void NoseNVTLeapfrogIntegrator::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("temperature",
               Value(myTemperature, ConstraintValueType::NotNegative()),
               Text("preferred system temperature")));
  parameters.push_back
    (Parameter("thermal", Value(myThermalInertia),
               Text( "heat bath coupling: 1.0 very, very strong, 0.0 none")));
  parameters.push_back
    (Parameter("bathPos", Value(myBathPosition), 0.0,
               Text("history of the difference of system and heat bath")));
}

STSIntegrator *NoseNVTLeapfrogIntegrator::doMake(const vector<Value> &values,
                                                 ForceGroup *fg) const {
  return new NoseNVTLeapfrogIntegrator(values[0], values[1], values[2],
                                       values[3], fg);
}
