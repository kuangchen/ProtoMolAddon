#include <protomol/integrator/normal/NormalModeLangLf.h>
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

using std::string;
using std::vector;

namespace ProtoMol {
  //__________________________________________________ NormalModeLangLf

  const string NormalModeLangLf::keyword( "NormalModeLangevinLeapfrog" );

  NormalModeLangLf::NormalModeLangLf() : MTSIntegrator(), NormalModeUtilities()
  {
  }

  NormalModeLangLf::NormalModeLangLf(int cycles, int firstmode, int nummode, Real gamma, int seed, Real temperature, bool gencn,
                     ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator) 
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), 
        NormalModeUtilities( firstmode, nummode, gamma, seed, temperature), genCompNoise(gencn)    
  {
  }


  NormalModeLangLf::~NormalModeLangLf() 
  {  
  }

  void NormalModeLangLf::initialize(ProtoMolApp *app){
    MTSIntegrator::initialize(app);
    //check valid eigenvectors
    //NM initialization if OK
    int nm_flags = NO_NM_FLAGS;
    if(genCompNoise) nm_flags |= GEN_COMP_NOISE;
    NormalModeUtilities::initialize((int)app->positions.size(), app, myForces, nm_flags); //last int for no complimentary forces or gen noise: GEN_COMP_NOISE
    //
    //do first force calculation, and remove non sub-space part
    app->energies.clear();	//Need this or initial error, due to inner integrator energy?
    initializeForces();
    //
    //take initial C velocites from system and remove non-subspace part
    if(*Q != NULL) subspaceVelocity(&app->velocities, &app->velocities);
    //

  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  void NormalModeLangLf::run(int numTimesteps) {
    //Real h = getTimestep() * Constant::INV_TIMEFACTOR;
    Real actTime;

    if( numTimesteps < 1 )
      return;

    //check valid eigenvectors
    if(*Q == NULL)
        report << error << "No Eigenvectors for NormalMode integrator."<<endr;
    //
    //time calculated in forces! so fix here
    actTime = app->topology->time + numTimesteps * getTimestep();
    //
    //main loop
    for( int i = 0; i < numTimesteps; i++ ) {
      //****main loop*************************************
      preStepModify();
      genProjGauss(&gaussRandCoord1, app->topology);
      doHalfKick();
      //
      //nmlDrift(&app->positions, &app->velocities, h, app->topology);
      doDrift();
      //constraints?
      app->energies.clear();
      //run minimizer if any remaining modes
      if(testRemainingModes()) myNextIntegrator->run(myCycleLength); //cyclelength 
      if(app->eigenInfo.reDiagonalize){	//rediagonalize?
            app->topology->time = actTime + (i - numTimesteps) * getTimestep();
            if(myPreviousIntegrator == NULL) 
                report << error << "[NormalModeLangLf::Run] Re-diagonalization forced with NormalModeLangLf as outermost Integrator. Aborting."<<endr;
            return;
      }
      //calculate sub space forces
      app->energies.clear();
      calculateForces();
      //
      genProjGauss(&gaussRandCoord1, app->topology);
      doHalfKick();
      //
      postStepModify();
    }
    //fix time
    app->topology->time = actTime;
    //
  }  

  void NormalModeLangLf::doDrift() {
      const Real h = getTimestep() * Constant::INV_TIMEFACTOR;
      app->positions.intoWeightedAdd(h, app->velocities);
      buildMolecularCenterOfMass(&app->positions, app->topology);
  }

  void NormalModeLangLf::doHalfKick() {
    const Real dt = getTimestep() * Constant::INV_TIMEFACTOR; // in fs
    const Real fdt = ( 1.0 - exp( -0.5 * myGamma * dt ) ) / myGamma;
    const Real vdt = exp(-0.5*myGamma*dt);
    const Real ndt = sqrt( ( 1.0 - exp( -myGamma * dt ) ) / (2.0 * myGamma) ); //was sqrt( fdt );
    const Real sqrtFCoverM = sqrt( 2.0 * Constant::BOLTZMANN * myTemp * myGamma );

    for( int i = 0; i < _N; i++ ) {
        // semi-update velocities
        app->velocities[i] = app->velocities[i]*vdt
                                +(*myForces)[i] * fdt / app->topology->atoms[i].scaledMass
                                    +gaussRandCoord1[i]*sqrtFCoverM*ndt;
    }
    subspaceVelocity(&app->velocities, &app->velocities);
    buildMolecularMomentum(&app->velocities, app->topology);
  }

  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeLangLf::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NotNegative()),1,Text("First mode to use in set")));
    parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NotNegative()),1,Text("Number of modes propagated")));
    parameters.push_back(Parameter("gamma",Value(myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
    parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234,Text("Langevin random seed")));
    parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
    parameters.push_back(Parameter("gencompnoise",Value(genCompNoise,ConstraintValueType::NoConstraints()),false,Text("Generate complimentary noise")));
 }

  MTSIntegrator* NormalModeLangLf::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeLangLf(values[0],values[1],values[2],values[3],values[4],values[5],values[6],fg,nextIntegrator);
  }

  void NormalModeLangLf::addModifierAfterInitialize(){
    adoptPostForceModifier(new ModifierForceProjection(this));
    MTSIntegrator::addModifierAfterInitialize();
  }

}
