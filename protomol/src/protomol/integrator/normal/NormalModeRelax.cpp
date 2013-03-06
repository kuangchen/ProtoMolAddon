#include <protomol/integrator/normal/NormalModeRelax.h>
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

using std::string;
using std::vector;

namespace ProtoMol {
  //__________________________________________________ NormalModeRelax

  const string NormalModeRelax::keyword( "NormalModeRelax" );

  NormalModeRelax::NormalModeRelax() : MTSIntegrator(), NormalModeUtilities()
  {
  }

  NormalModeRelax::NormalModeRelax(int cycles, Real minimlim, bool rediag, bool simplemin, 
                                    ForceGroup *overloadedForces, StandardIntegrator *nextIntegrator)
    : MTSIntegrator(cycles, overloadedForces, nextIntegrator), 
        NormalModeUtilities( 1, 1, 91.0, 1234, 300.0),
            minLim(minimlim), reDiag(rediag), simpleMin(simplemin)
  {
  }


  NormalModeRelax::~NormalModeRelax() 
  {  
  }

  void NormalModeRelax::initialize(ProtoMolApp* app){
    MTSIntegrator::initialize(app);
    //test not topmost integrator
    if(top() == this) report << error << "NormalModeRelax cannot be top integrator."<<endr;
    //
    initializeForces();
    //
    myPreviousNormalMode  = dynamic_cast<NormalModeUtilities*>(myPreviousIntegrator);
    //check valid eigenvectors
    firstMode = myPreviousNormalMode->firstMode; numMode = myPreviousNormalMode->numMode;
    //NM initialization if OK
    NormalModeUtilities::initialize((int)app->positions.size(), app, myForces, COMPLIMENT_FORCES); //last for complimentary forces
    //
    app->energies.clear();	//Need this or initial error, due to inner integrator energy?
    //
    //diagnostics
    avItrs = 0;		//average number of minimizer iterations/force calcs
    avMinForceCalc = 0;
    numSteps = 0;	//total steps

  }

  //*************************************************************************************
  //****Normal run routine***************************************************************
  //*************************************************************************************

  void NormalModeRelax::run(int numTimesteps) {

    if( numTimesteps < 1 )
      return;

    //check valid eigenvectors
    if(*Q == NULL)
        report << error << "No Eigenvectors for NormalMode integrator."<<endr;
    //
    //do minimization with local forces, max loop 100, set subSpace minimization true
    itrs = minimizer(minLim, 100, simpleMin, reDiag, true, &forceCalc, &lastLambda, &app->energies, &app->positions, app->topology);
    Real minPotEnergy = app->energies.potentialEnergy();				//save potential energy before random
    //constraints?
    app->energies.clear();
    //run brownian if any remaining modes
    if(testRemainingModes()) myNextIntegrator->run(myCycleLength);	//cyclelength 
    (app->energies)[ScalarStructure::OTHER] = minPotEnergy;			//store minimum potential energy
    if(*Q == NULL){	//rediagonalize?
        if(myPreviousIntegrator == NULL) 
            report << error << "[NormalModeRelax::Run] Re-diagonalization forced with NormalModeRelax as outermost Integrator. Aborting."<<endr;
        return;
    }
    //
    postStepModify();
    //
  }  

  //*************************************************************************************
  //****Output int paramiters************************************************************
  //*************************************************************************************

  void NormalModeRelax::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("minimlim",Value(minLim,ConstraintValueType::NotNegative()),0.1,Text("Minimizer target PE difference kcal mole^{-1}")));
    parameters.push_back(Parameter("rediag",Value(reDiag,ConstraintValueType::NoConstraints()),false,Text("Force re-digonalize")));
    parameters.push_back(Parameter("simplemin",Value(simpleMin,ConstraintValueType::NoConstraints()),true,Text("Simple minimizer or exact minima projection.")));
  }

  MTSIntegrator* NormalModeRelax::doMake(const vector<Value>& values,ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    return new NormalModeRelax(values[0],values[1],values[2],values[3],fg,nextIntegrator);
  }

  //void NormalModeRelax::addModifierAfterInitialize(){
  //  adoptPostForceModifier(new ModifierForceProjection(this));
  //  MTSIntegrator::addModifierAfterInitialize();
  //}

  //*************************************************************************************
  //****Minimizers virtual force calculation*********************************************
  //*************************************************************************************

  void NormalModeRelax::utilityCalculateForces(){
      app->energies.clear();	//need this as MTS, not innermost
      calculateForces();
  }

}
