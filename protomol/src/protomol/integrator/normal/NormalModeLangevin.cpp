#include <protomol/integrator/normal/NormalModeLangevin.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/integrator/normal/ModifierForceProjection.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ NormalModeLangevin

const string NormalModeLangevin::keyword("NormalModeLangevin");

NormalModeLangevin::NormalModeLangevin() :
  MTSIntegrator(), NormalModeUtilities()
{}

NormalModeLangevin::NormalModeLangevin(int cycles, int firstmode, int nummode,
                                       Real gamma, int seed, Real temperature,
                                       bool gencn,
                                       ForceGroup *overloadedForces,
                                       StandardIntegrator *nextIntegrator) :
  MTSIntegrator(cycles, overloadedForces,
                nextIntegrator),
  NormalModeUtilities(firstmode, nummode, gamma, seed,
                      temperature),
  genCompNoise(gencn)
{}

NormalModeLangevin::~NormalModeLangevin()
{}

void NormalModeLangevin::initialize(ProtoMolApp *app) {
  MTSIntegrator::initialize(app);
  //check valid eigenvectors
  //NM initialization if OK
  int nm_flags = NO_NM_FLAGS;
  if (genCompNoise) nm_flags |= GEN_COMP_NOISE;
  
  //last int for no complimentary forces or gen noise: GEN_COMP_NOISE
  NormalModeUtilities::initialize((int)app->positions.size(), app,
                                  myForces, nm_flags);
  //
  //do first force calculation, and remove non sub-space part
  //Need this or initial error, due to inner integrator energy?
  app->energies.clear();
  initializeForces();
  //
  //take initial C velocites from system and remove non-subspace part
  if (*Q != 0) subspaceVelocity(&app->velocities, &app->velocities);
  //
}

//******************************************************************************
//****Normal run routine********************************************************
//******************************************************************************

void NormalModeLangevin::run(int numTimesteps) {
  Real h = getTimestep() * Constant::INV_TIMEFACTOR;
  Real actTime;

  if (numTimesteps < 1)
    return;

  //check valid eigenvectors
  if (*Q == 0)
    report << error << "No Eigenvectors for NormalMode integrator." << endr;
  //
  //time calculated in forces! so fix here
  actTime = app->topology->time + numTimesteps *getTimestep();

  //
  //main loop
  for (int i = 0; i < numTimesteps; i++) {
    //****main loop*************************************
    preStepModify();
    doHalfKick();
    //
    nmlDrift(&app->positions, &app->velocities, h, app->topology);
    //constraints?
    app->energies.clear();
    //run minimizer if any remaining modes
    //cyclelength
    if (testRemainingModes()) myNextIntegrator->run(myCycleLength);
    if (app->eigenInfo.reDiagonalize) {   //rediagonalize?
      app->topology->time = actTime + (i - numTimesteps) * getTimestep();
      if (myPreviousIntegrator == 0)
        report
          << error << "[NormalModeLangevin::Run] Re-diagonalization forced "
          "with NormalModeLangevin as outermost Integrator. Aborting."
          << endr;
      return;
    }
    //calculate sub space forces
    app->energies.clear();
    calculateForces();
    //
    doHalfKick();
    //
    postStepModify();
  }

  //fix time
  app->topology->time = actTime;
  //
}

//******************************************************************************
//****Output int paramiters*****************************************************
//******************************************************************************

void NormalModeLangevin::getParameters(vector<Parameter> &parameters) const {
  MTSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("firstmode",
               Value(firstMode, ConstraintValueType::NotNegative()), 1,
               Text("First mode to use in set")));
  parameters.push_back
    (Parameter("numbermodes",
               Value(numMode, ConstraintValueType::NotNegative()), 1,
               Text("Number of modes propagated")));
  parameters.push_back
    (Parameter("gamma", Value(myGamma * (1000 * Constant::INV_TIMEFACTOR),
                              ConstraintValueType::NotNegative()),
               80.0, Text("Langevin Gamma")));
  parameters.push_back
    (Parameter("seed", Value(mySeed, ConstraintValueType::NotNegative()),
               1234, Text("Langevin random seed")));
  parameters.push_back
    (Parameter("temperature", Value(myTemp, ConstraintValueType::NotNegative()),
               300.0, Text("Langevin temperature")));
  parameters.push_back
    (Parameter("gencompnoise",
 Value(genCompNoise, ConstraintValueType::NoConstraints()),
               false, Text("Generate complimentary noise")));
}

MTSIntegrator *NormalModeLangevin::doMake(const vector<Value> &values,
                                          ForceGroup *fg,
                                          StandardIntegrator *nextIntegrator)
const {
  return new NormalModeLangevin(values[0], values[1], values[2], values[3],
                                values[4], values[5], values[6], fg,
                                nextIntegrator);
}

void NormalModeLangevin::addModifierAfterInitialize() {
  adoptPostForceModifier(new ModifierForceProjection(this));
  MTSIntegrator::addModifierAfterInitialize();
}
