#include <protomol/integrator/normal/NormalModeMinimizer.h>
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>
#include <protomol/ProtoMolApp.h>

//#include "ModifierForceProjection.h"

using namespace std;

using namespace ProtoMol::Report;

using std::string;
using std::vector;


namespace ProtoMol {
  //__________________________________________________ NormalModeMinimizer

  const string NormalModeMinimizer::keyword( "NormalModeMinimizer" );

  NormalModeMinimizer::NormalModeMinimizer() :
    avItrs(0), itrs(0), avMinForceCalc(0), numSteps(0), minCount(0),
    forceCalc(0), minLim(0), randforce(0), myPreviousNormalMode(0),
    lastLambda(0), reDiag(0), simpleMin(0), eUFactor(0), randStp(0),
    rediagOnMaxMinSteps(0) {}

  NormalModeMinimizer::NormalModeMinimizer
  (Real timestep, int firstmode, int nummode, Real gamma, int seed,
   Real temperature,  Real minimlim, int rforce, bool red, bool simplemin,
   int redmaxmin, ForceGroup *overloadedForces) 
    : STSIntegrator(timestep, overloadedForces),
      NormalModeUtilities(firstmode, nummode, gamma, seed, temperature),
      avItrs(0), itrs(0), avMinForceCalc(0), numSteps(0), minCount(0),
      forceCalc(0), minLim(minimlim), randforce(rforce),
      myPreviousNormalMode(0), lastLambda(0), reDiag(red), 
      simpleMin(simplemin), eUFactor(0), randStp(0),
      rediagOnMaxMinSteps(redmaxmin) {}

  NormalModeMinimizer::~NormalModeMinimizer() 
  {  
      report.precision(5);
      if(numSteps) report <<plain<<"NML Stats.: Minimizer average iterations = "<<(float)avItrs/(float)numSteps
          <<", force average calcs = "<<(float)avMinForceCalc/(float)numSteps<<endl;

  }

  void NormalModeMinimizer::initialize(ProtoMolApp* app){
    STSIntegrator::initialize(app);
    initializeForces();
    myPreviousNormalMode  = dynamic_cast<NormalModeUtilities*>(myPreviousIntegrator);
    //Check if using complement of prev integrator, then copy all integrator paramiters
    if(firstMode == -1 && numMode == -1){
        firstMode = myPreviousNormalMode->firstMode; numMode = myPreviousNormalMode->numMode;
        myGamma = myPreviousNormalMode->myGamma; mySeed = myPreviousNormalMode->mySeed; myTemp = myPreviousNormalMode->myTemp;
    }
    //NM initialization if OK
    NormalModeUtilities::initialize((int)app->positions.size(), app,
				    myForces, COMPLIMENT_FORCES); //last for complimentary forces, no gen noise
    //
    //initialize minimizer noise vars
    randStp = 0.0;
    //***********************************************************	
    //diagnostics
    avItrs = 0;		//average number of minimizer iterations/force calcs
    avMinForceCalc = 0;
    numSteps = 0;	//total steps

    //Set up minimum limit
    app->eigenInfo.myMinimumLimit = minLim;
    //
    gaussRandCoord1.zero(-1); //zero vector for use as random force
  }

  void NormalModeMinimizer::run(int numTimesteps) {
    if( numTimesteps < 1 )
        return;

    //check valid eigenvectors
    if(*Q == NULL)
        report << error << "No Eigenvectors for NormalMode integrator."<<endr;
    //
    preStepModify();
    //remove last random pertubation
    if(randforce) app->positions.intoSubtract(gaussRandCoord1);
    //(*myPositions).intoWeightedAdd(-randStp,gaussRandCoord1);
    //do minimization with local forces, max loop rediagOnMaxMinSteps, set subSpace minimization true
    itrs = minimizer(minLim, rediagOnMaxMinSteps>0 ? rediagOnMaxMinSteps : 100 , simpleMin, reDiag, true, &forceCalc, &lastLambda, &app->energies, &app->positions, app->topology);

    //flag excessive minimizations
    if(itrs > 10) report << debug(1) << "[NormalModeMinimizer::run] iterations = " << itrs << "." << endr;
    else report << debug(3) << "[NormalModeMinimizer::run] iterations = " << itrs << "." << endr;

    avItrs += itrs;

    //rediagonalize if minimization steps exceeds 'rediagOnMaxMinSteps'
    if(reDiag && rediagOnMaxMinSteps > 0 && itrs >= rediagOnMaxMinSteps){
      report << debug(1) << "[NormalModeMinimizer::run] Minimization steps (" 
              << itrs << ") exceeded maximum (" << rediagOnMaxMinSteps << "), forcing re-diagonalize." << endr;
      itrs = -1;  //force re-diag
    }

    numSteps++;
    avMinForceCalc += forceCalc;
    report <<debug(5)<<"[NormalModeMinimizer::run] iterations = "<<itrs<<" average = "<<
                (float)avItrs/(float)numSteps<<" force calcs = "<<forceCalc<<" average = "<<(float)avMinForceCalc/(float)numSteps<<endl;
    if(randforce && itrs != -1 && lastLambda > 0){	//add random force, but not if rediagonalizing

      //add random force
      if(myPreviousNormalMode->genCompNoise) gaussRandCoord1 = myPreviousNormalMode->tempV3DBlk;
      else genProjGauss(&gaussRandCoord1, app->topology);
      //lambda = eUFactor / *eigValP;	//user factor
      randStp = sqrt(2 * Constant::BOLTZMANN * myTemp * lastLambda);	//step length
      //(*myPositions).intoWeightedAdd(randStp,gaussRandCoord1);
      gaussRandCoord1.intoWeighted(randStp,gaussRandCoord1);
      app->positions.intoAdd(gaussRandCoord1);

      //additional random steps?
      if(randforce > 1){
        eUFactor = 0.5;
        Real lambda;
        for(int i=1;i<randforce;i++){
          utilityCalculateForces();
          nonSubspaceForce(myForces, myForces);
          for(int j=0;j<_N;j++) (*myForces)[j] /= app->topology->atoms[j].scaledMass;
          lambda = eUFactor / *eigValP;	//user factor
          //update positions
          app->positions.intoWeightedAdd(lambda,*myForces);
          gaussRandCoord1.intoWeightedAdd(lambda,*myForces); //add to grc1 so can be removed at the next step
          //random force
          genProjGauss(&gaussRandCoord2, app->topology);
          randStp = sqrt(2 * Constant::BOLTZMANN * myTemp * lambda);
          app->positions.intoWeightedAdd(randStp,gaussRandCoord2);
          //add to grc1 so can be removed at the next step
          gaussRandCoord1.intoWeightedAdd(randStp,gaussRandCoord2); 
        }
      }

    }else{
        gaussRandCoord1.zero(-1);

        //flag re-diagonalize if detected
        if(itrs == -1) app->eigenInfo.reDiagonalize = true;
    }

    //
    postStepModify();
  }  

  void NormalModeMinimizer::getParameters(vector<Parameter>& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("firstmode",Value(firstMode,ConstraintValueType::NoConstraints()),-1,Text("First mode to use in set")));
    parameters.push_back(Parameter("numbermodes",Value(numMode,ConstraintValueType::NoConstraints()),-1,Text("Number of modes propagated")));
    parameters.push_back(Parameter("gamma",Value(myGamma*(1000 * Constant::INV_TIMEFACTOR),ConstraintValueType::NotNegative()),80.0,Text("Langevin Gamma")));
    parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234,Text("Langevin random seed")));
    parameters.push_back(Parameter("temperature",Value(myTemp,ConstraintValueType::NotNegative()),300.0,Text("Langevin temperature")));
    parameters.push_back(Parameter("minimlim",Value(minLim,ConstraintValueType::NotNegative()),0.1,Text("Minimizer target PE difference kcal mole^{-1}")));
    parameters.push_back(Parameter("randforce",Value(randforce,ConstraintValueType::NotNegative()),1,Text("Add random force/EM steps")));
    parameters.push_back(Parameter("rediag",Value(reDiag,ConstraintValueType::NoConstraints()),false,Text("Force re-digonalize")));
    parameters.push_back(Parameter("simplemin",Value(simpleMin,ConstraintValueType::NoConstraints()),true,Text("Simple minimizer or exact minima projection.")));
    parameters.push_back(Parameter("rediagmaxminsteps",Value(rediagOnMaxMinSteps,ConstraintValueType::NotNegative()),0,Text("Rediagonalize if maximum minimizer steps exceeded.")));

  }

  STSIntegrator* NormalModeMinimizer::doMake(const vector<Value>& values,ForceGroup* fg)const{
    return new NormalModeMinimizer(values[0],values[1],values[2],values[3],values[4],values[5],
                                   values[6],values[7],values[8],values[9],values[10],fg);
  }

  //void NormalModeMinimizer::addModifierAfterInitialize(){
  //  adoptPostForceModifier(new ModifierForceProjection(this));
  //  STSIntegrator::addModifierAfterInitialize();
  //}

  //*************************************************************************************
  //****Minimizers virtual force calculation*********************************************
  //*************************************************************************************

  void NormalModeMinimizer::utilityCalculateForces(){
      calculateForces();
  }

}

