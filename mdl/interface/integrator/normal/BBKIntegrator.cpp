#include "BBKIntegrator.h"
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/base/PMConstants.h>


using namespace ProtoMol::Report;
using std::vector;
using std::string;

namespace ProtoMol {
  //________________________________________________________ BBKIntegrator

  const string BBKIntegrator::keyword("BBK");


  BBKIntegrator::BBKIntegrator():
    STSIntegrator(),
    myLangevinTemperature(-1.0),
    myGamma(-1.0), 
    mySeed(-1){}

  BBKIntegrator::BBKIntegrator(Real timestep,
			       Real langevinTemperature,
			       Real gamma,
			       int seed,
			       ForceGroup *overloadedForces):
    STSIntegrator(timestep,overloadedForces),
    myLangevinTemperature(langevinTemperature),
    myGamma(gamma*0.001), 
    mySeed(seed){
  }

  void BBKIntegrator::initialize(ProtoMolApp* app){
  
    STSIntegrator::initialize(app);
    initializeForces();
  }

  void BBKIntegrator::run (int numTimesteps) {

    for (int i = 0; i < numTimesteps; i++) {
      preStepModify();
      addRandomForceToMyForce();
      doFirstHalfKick();
      doDriftOrNextIntegrator();
      calculateForces();
      addRandomForceToMyForce();
      doSecondHalfKick();
      postStepModify();
    }  
  }



  // BBK procedure to add random force to the force vector
  void BBKIntegrator::addRandomForceToMyForce() {
    const Real dt = getTimestep()/Constant::TIMEFACTOR;
    Real forceConstant =(2*Constant::BOLTZMANN*myLangevinTemperature*myGamma)/dt;
    for(unsigned int i=0;i<app->positions.size();i++){
      //  Generate gaussian random numbers for each spatial directions such that 
      //  the average is 0 and the standard deviation is fParam.  The algorithm 
      //  that is  used was taken from X-PLOR. The algorithm is a "sum of 
      //  uniform deviates algorithm" which may be found in Abramowitz and 
      //  Stegun, "Handbook of Mathematical Functions", pg 952.
      Vector3D randomForce(randomGaussianNumber(mySeed),randomGaussianNumber(mySeed),randomGaussianNumber(mySeed));
      Real mass = app->topology->atoms[i].scaledMass;
      randomForce *= sqrt(forceConstant*mass);
      (*myForces)[i] += randomForce;
    }
  }

  void BBKIntegrator::doFirstHalfKick(){
    // Our range. This also assumed by the methods in Parallel!!!
    Real dt = getTimestep()/Constant::TIMEFACTOR;
    // We only compute the doKick & doDrift on our range
    for(unsigned int i=0;i<app->positions.size();i++){
      (app->velocities)[i] *= 1.0-0.5*dt*myGamma;
      (app->velocities)[i] += (*myForces)[i]*dt*0.5/app->topology->atoms[i].scaledMass;
    }
    buildMolecularMomentum(&app->velocities,app->topology);
  }

  void BBKIntegrator::doSecondHalfKick(){
    // Our range. This also assumed by the methods in Parallel!!!
    Real dt = getTimestep()/Constant::TIMEFACTOR;
    // We only compute the doKick & doDrift on our range
    for(unsigned int i=0;i<app->positions.size();i++){
      (app->velocities)[i] += (*myForces)[i]*dt*0.5/app->topology->atoms[i].scaledMass;
      (app->velocities)[i] /= (1.0+0.5*dt*myGamma);
    }
    buildMolecularMomentum(&app->velocities,app->topology);
  }

  void BBKIntegrator::getParameters(vector<Parameter>& parameters) const {
    STSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("temperature",Value(myLangevinTemperature,ConstraintValueType::NotNegative())));
    parameters.push_back(Parameter("gamma",Value(myGamma*1000,ConstraintValueType::NotNegative())));
    parameters.push_back(Parameter("seed",Value(mySeed,ConstraintValueType::NotNegative()),1234));
  }

  STSIntegrator* BBKIntegrator::doMake(const vector<Value>& values,ForceGroup* fg)const{
    return new BBKIntegrator(values[0],values[1],values[2],values[3],fg);

  }

}
