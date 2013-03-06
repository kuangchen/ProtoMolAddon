#include <protomol/integrator/normal/NormalModeMori.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

#include <protomol/integrator/normal/ModifierForceProjection.h>

#include <fstream>
using namespace std;


using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol {
  //__________________________________________________ NormalModeMori

  const string NormalModeMori::keyword( "NormalModeMori" );

  NormalModeMori::NormalModeMori() : MTSIntegrator(), NormalModeUtilities()
  {
      ex0=NULL;//####diagnostics
  }

  NormalModeMori::NormalModeMori(int cycles, int firstmode, int nummode, Real gamma, int seed, Real temperature,
                    std::string mof, bool instf, //####added mof for diagnostics, instf for NML comparison
                        ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator) 
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), 
        NormalModeUtilities( firstmode, nummode, gamma, seed, temperature),     
            modeOutput(mof), instForce(instf)//####added modeOutput for diagnostics, instForce for NML comparison

  {
      ex0=NULL;//####diagnostics
  }


  NormalModeMori::~NormalModeMori() 
  {        
      if(ex0!=NULL) delete ex0;//####diagnostics
  }

  void NormalModeMori::initialize(ProtoMolApp* app){
    MTSIntegrator::initialize(app);
    //
    //point to bottom integrator, for average force
    myBottomNormalMode  = dynamic_cast<NormalModeUtilities*>(bottom());
    //check valid eigenvectors
    //NM initialization if OK
    NormalModeUtilities::initialize((int)app->positions.size(), app, myForces, NO_NM_FLAGS); //last for non-complimentary forces
    //
    //do first force calculation, and remove non sub-space part
    app->energies.clear();	//Need this or initial error, due to inner integrator energy?
    initializeForces();
    //
    //take initial C velocites from system and remove non-subspace part
    if(*Q != NULL) subspaceVelocity(&app->velocities, &app->velocities);
    //	
    //####diagnostics
    if(modeOutput != ""){
        //save initial positions
        ex0 = new Vector3DBlock;
        *ex0 = app->positions;
        ofstream myFile;
        myFile.open(modeOutput.c_str(),ofstream::out);
        //close file
        myFile.close();
    }

  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  void NormalModeMori::run(int numTimesteps) {
    Real h = getTimestep() * Constant::INV_TIMEFACTOR;
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
      doHalfKick();
      //
      nmlDrift(&app->positions, &app->velocities, h, app->topology);
      //constraints?
      app->energies.clear();
      //run minimizer if any remaining modes
      if(testRemainingModes()) myNextIntegrator->run(myCycleLength); //cyclelength 
      if(*Q == NULL){	//rediagonalize?
            app->topology->time = actTime - (i - numTimesteps) * getTimestep();
            if(myPreviousIntegrator == NULL) 
                report << error << "[NormalModeMori::Run] Re-diagonalization forced with NormalModeMori as outermost Integrator. Aborting."<<endr;
            return;
      }
      //#################Put averaged force code here##############################
      //calculate sub space forces, just do this for the energy
      Real minPotEnergy = (app->energies)[ScalarStructure::OTHER];	//save minimum potential energy
      app->energies.clear();
      calculateForces();
      //transfer the real average force, Note: force fields must be identical
      if(!instForce){
        for( unsigned int i=0;i < app->positions.size(); i++) (*myForces)[i] = myBottomNormalMode->tempV3DBlk[i];
      }
      (app->energies)[ScalarStructure::OTHER] = minPotEnergy;			//restore minimum potential energy
      //###########################################################################
      //
      doHalfKick();
      //
      postStepModify();
    }	
    //####diagnostics
    //output modes?
    if(modeOutput != ""){
        //get modes
        modeSpaceProj(tmpC, &app->positions, ex0);
        //Output modes for analysis
        ofstream myFile;
        myFile.open(modeOutput.c_str(),ofstream::app);
        myFile.precision(10);
        for(int ii=0;ii<_rfM-1;ii++) myFile << tmpC[ii] << ", ";
        myFile << tmpC[_rfM-1] << endl; //last
        //close file
        myFile.close();
    }
    //####
    //fix time
    app->topology->time = actTime;
    //
  }  

  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeMori::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NotNegative()),1,Text("First mode to use in set")));
    parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NotNegative()),1,Text("Number of modes propagated")));
    parameters.push_back(Parameter("gamma",Value(myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
    parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234,Text("Langevin random seed")));
    parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
    //####diagnostics
    parameters.push_back(Parameter("modeOutput",Value(modeOutput,ConstraintValueType::NoConstraints()),std::string(""),Text("Mode output filename")));
    parameters.push_back(Parameter("instForce",Value(instForce,ConstraintValueType::NoConstraints()),false,Text("Use instantaneous force")));
    //####
  }

  MTSIntegrator* NormalModeMori::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeMori(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],fg,nextIntegrator);
  }

  void NormalModeMori::addModifierAfterInitialize(){
    adoptPostForceModifier(new ModifierForceProjection(this));
    MTSIntegrator::addModifierAfterInitialize();
  }

}
