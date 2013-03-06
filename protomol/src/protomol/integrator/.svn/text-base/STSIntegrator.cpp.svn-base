#include <protomol/integrator/STSIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/modifier/ModifierIncrementTimestep.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>

using namespace ProtoMol::Report;
using namespace std;
using namespace ProtoMol;
//____ STSIntegrator

STSIntegrator::STSIntegrator() :
  StandardIntegrator(), myTimestep(0.0) {}

STSIntegrator::STSIntegrator(Real timestep, ForceGroup *overloadedForces) :
  StandardIntegrator(overloadedForces), myTimestep(timestep) {}

void STSIntegrator::addModifierAfterInitialize() {
  adoptPostForceModifier(new ModifierIncrementTimestep());
}

void STSIntegrator::doDriftOrNextIntegrator() {
  preDriftOrNextModify();
  doDrift();
  postDriftOrNextModify();
}

void STSIntegrator::calculateForces() {
  //  -------------------------------------------------------------------  //
  //  Energies have to be cleared for the innermost integrator which is    //
  //  always an STS integrator.  Potential energy is U^{short}(X1) +       //
  //  U^{long}(X1) + ... +U^{last}(X1) forces are gradient of potential,   //
  //  thus they are set to zero here also.                                 //
  //                                                                       //
  //  Forces are cleared in StandardIntegrator::calculateForces()          //
  //  -------------------------------------------------------------------  //

  app->energies.clear();
  StandardIntegrator::calculateForces();
}

void STSIntegrator::doDrift() {
  Real h = getTimestep() * Constant::INV_TIMEFACTOR;
  app->positions.intoWeightedAdd(h, app->velocities);
  buildMolecularCenterOfMass(&app->positions, app->topology);
}

void STSIntegrator::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter("timestep", Value(myTimestep, ConstraintValueType::Positive())));
}

STSIntegrator *STSIntegrator::make(const vector<Value> &values,
                                   ForceGroup *fg) const {
  assertParameters(values);

  return adjustAlias(doMake(values, fg));
}
