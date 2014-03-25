#include <protomol/addon/damping/LeapfrogDampedIntegrator.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/addon/Constants.h>

#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <iterator>


using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
using std::runtime_error;
using std::istringstream;
using std::copy;
using std::inserter;
using namespace ProtoMolAddon::Damping;
namespace pc = ProtoMol::Constant;

//____ LeapfrogDampedIntegrator

const string LeapfrogDampedIntegrator::keyword("LeapfrogDamped");

LeapfrogDampedIntegrator::LeapfrogDampedIntegrator() :
  STSIntegrator(), spec_fname(), spec() {}

LeapfrogDampedIntegrator::LeapfrogDampedIntegrator(Real timestep,
						   const string &spec_fname,
						   ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces), spec_fname(spec_fname) {

  ifstream is(spec_fname);
  if (!is) 
    throw runtime_error("Error opening file " + spec_fname);

  is >> spec;
}

LeapfrogDampedIntegrator::~LeapfrogDampedIntegrator() {}

void LeapfrogDampedIntegrator::initialize(ProtoMolApp *app) {
    STSIntegrator::initialize(app);
    initializeForces();
}

void LeapfrogDampedIntegrator::doHalfKickdoDrift() {
    if (anyPreDriftOrNextModify()) {
	doHalfKick();
	doDriftOrNextIntegrator();
    } else {
	Real h = getTimestep() * pc::INV_TIMEFACTOR;
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

void LeapfrogDampedIntegrator::doKickdoDrift() {
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
      Real h = getTimestep() * pc::INV_TIMEFACTOR;
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

void LeapfrogDampedIntegrator::run(int numTimesteps) {
    if (numTimesteps < 1)
	return;
    preStepModify();
    doHalfKickdoDrift();
    calculateForces();
    for (int i = 1; i < numTimesteps; i++) {
	doKickdoDrift();
	calculateForces();
    }

    doHalfKick();
    postStepModify();
    
    Real h = getTimestep() * pc::INV_TIMEFACTOR;
    Real now = app->topology->time * ProtoMolAddon::Constant::TIME_CONV;
    spec.Damp(app, now, h);
}


STSIntegrator *LeapfrogDampedIntegrator::doMake(const vector<Value> &values,
						ForceGroup *fg) const {
  return new LeapfrogDampedIntegrator(values[0], values[1], fg);
}

//  --------------------------------------------------------------------  //
//  This function is necessary to compute the shadow Hamiltonian and it   //
//  is integrator specific.  This version is written to work with LF.     //
//  Update beta: beta -= dt * ( q * F + 2 U )                             //
//  --------------------------------------------------------------------  //

void LeapfrogDampedIntegrator::updateBeta(Real dt) {
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


void LeapfrogDampedIntegrator::getParameters(vector<Parameter> &parameters) const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back(Parameter("spec_fname", Value(spec_fname, ConstraintValueType::NotEmpty())));
}

