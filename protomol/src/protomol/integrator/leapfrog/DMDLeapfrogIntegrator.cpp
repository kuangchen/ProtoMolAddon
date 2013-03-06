#include <protomol/integrator/leapfrog/DMDLeapfrogIntegrator.h>
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
using namespace ProtoMol;

//____ DMDLeapfrogIntegrator

const string DMDLeapfrogIntegrator::keyword("DMDLeapfrog");

DMDLeapfrogIntegrator::DMDLeapfrogIntegrator() :
  STSIntegrator(), myDissipativeForces(0), myRandomForces(0), myVhat(0),
  myGamma(-1.0), myTemperature(-1.0), myNumIter(-1), mySigma(0.0), mySeed(-1)
{}

DMDLeapfrogIntegrator::DMDLeapfrogIntegrator(Real timestep, int numIter,
                                             Real gamma, Real temperature,
                                             int seed,
                                             ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces),
  myDissipativeForces(new Vector3DBlock),
  myRandomForces(new Vector3DBlock), myVhat(new Vector3DBlock),
  myGamma(gamma * 0.001), myTemperature(temperature), myNumIter(numIter),
  mySigma(sqrt(2 * myGamma * myTemperature * Constant::BOLTZMANN)), mySeed(seed)
{}

DMDLeapfrogIntegrator::~DMDLeapfrogIntegrator() {
  if (myDissipativeForces != 0) delete myDissipativeForces;
  if (myRandomForces != 0) delete myRandomForces;
  if (myVhat != 0) delete myVhat;
}

void DMDLeapfrogIntegrator::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();

  myDissipativeForces->zero(app->positions.size());
  myRandomForces->zero(app->positions.size());
  myVhat->zero(app->positions.size());
  calculateDissipativeAndRandomForces();
}

void DMDLeapfrogIntegrator::doHalfKick() {
  Real f = 0.5 * getTimestep() / Constant::TIMEFACTOR;
  Real g = 0.5 * sqrt(2 * f);
  Vector3DBlock myTempCoordBlock(app->positions.size());
  myTempCoordBlock.zero();
  myTempCoordBlock.intoWeightedAdd(f, *myForces);
  myTempCoordBlock.intoWeightedAdd(f, *myDissipativeForces);
  myTempCoordBlock.intoWeightedAdd(g, *myRandomForces);
  // We only compute the doHalfKick and doDrift in our range
  for (unsigned int i = 0; i < app->positions.size(); i++)
    app->velocities[i] +=
      myTempCoordBlock[i] / app->topology->atoms[i].scaledMass;

  buildMolecularMomentum(&app->velocities, app->topology);
}

void DMDLeapfrogIntegrator::doHalfKickVhat() {
  Real f = 0.5 * getTimestep() / Constant::TIMEFACTOR;
  Real g = 0.5 * sqrt(2 * f);
  Vector3DBlock myTempCoordBlock(app->positions.size());
  myTempCoordBlock.zero();
  myTempCoordBlock.intoWeightedAdd(f, *myForces);
  myTempCoordBlock.intoWeightedAdd(g, *myRandomForces);
  // We only compute the doHalfKick and doDrift in our range
  for (unsigned int i = 0; i < app->positions.size(); i++)
    (*myVhat)[i] = app->velocities[i] +
                   myTempCoordBlock[i] / app->topology->atoms[i].scaledMass;
}

void DMDLeapfrogIntegrator::doHalfKickIterate() {
  Real f = 0.5 * getTimestep() / Constant::TIMEFACTOR;
  Vector3DBlock myTempCoordBlock(app->positions.size());
  Vector3DBlock vHat(app->positions.size());
  vHat.zero();
  vHat.intoAdd(app->velocities);
  for (int j = 0; j < myNumIter; j++) {
    myTempCoordBlock.zero();
    myTempCoordBlock.intoWeightedAdd(f, *myDissipativeForces);
    // We only compute the doHalfKick and doDrift in our range
    for (unsigned int i = 0; i < app->positions.size(); i++)
      app->velocities[i] = 
        (*myVhat)[i] + myTempCoordBlock[i] / app->topology->atoms[i].scaledMass;

    calculateDissipativeForces();
  }

  buildMolecularMomentum(&app->velocities, app->topology);
}

//Here I did a little different than the litterature:
//"Towards better integrators for dissipative particle dynamics simulations"
//by Gerhard Besold, etc. In their paper, they defined rij=ri-rj and vij=vi-vj,
//which does not conform with mathematical convention. I changed them into
//rij=rj-ri, and vij=vj-vi. Thus the dissipitive pairwise force F_ij,
//F_ij = - gamma*weight^2*(v_ij \dot e_ij) e_ij, where e_ij is the unit vector.
//This is the force exerted on particle "j" by particle "i".
void DMDLeapfrogIntegrator::calculateDissipativeForces() {
  myDissipativeForces->zero();
  for (unsigned int i = 0; i < app->topology->angles.size(); i++) {
    int a1 = app->topology->angles[i].atom1;
    int a2 = app->topology->angles[i].atom2;
    int a3 = app->topology->angles[i].atom3;

    const Vector3D &posi1 = app->positions[a1];
    const Vector3D &posi2 = app->positions[a2];
    const Vector3D &posi3 = app->positions[a3];
    const Vector3D &vel1 = app->velocities[a1];
    const Vector3D &vel2 = app->velocities[a2];
    const Vector3D &vel3 = app->velocities[a3];

    //now handle atoms 1 and 2
    Vector3D unitVec = posi2 - posi1;   // now it is only a vector
    unitVec.normalize();
    // now the vector is unit vector, the dist is the length of the original 
    // vector.  We assume the weight, w1, and w2 to be 1 always, thus the
    // following still holds true
    // w1 = w2^2;
    Real coeff = -myGamma * ((vel2 - vel1).dot(unitVec));
    (*myDissipativeForces)[a2] += unitVec * coeff;
    (*myDissipativeForces)[a1] -= unitVec * coeff;

    //now handle atoms 2 and 3
    unitVec = posi3 - posi2;   // now it is only a vector
    unitVec.normalize();

    coeff = -myGamma * ((vel3 - vel2).dot(unitVec));
    (*myDissipativeForces)[a3] += unitVec * coeff;
    (*myDissipativeForces)[a2] -= unitVec * coeff;

    //now handle atoms 1 and 3
    unitVec = posi3 - posi1;   // now it is only a vector
    unitVec.normalize();
    coeff = -myGamma * ((vel3 - vel1).dot(unitVec));
    (*myDissipativeForces)[a3] += unitVec * coeff;
    (*myDissipativeForces)[a1] -= unitVec * coeff;
  }
}

void DMDLeapfrogIntegrator::calculateDissipativeAndRandomForces() {
  Real sdv = 1;
  myRandomForces->zero();
  myDissipativeForces->zero();
  for (unsigned int i = 0; i < app->topology->angles.size(); i++) {
    int a1 = app->topology->angles[i].atom1;
    int a2 = app->topology->angles[i].atom2;
    int a3 = app->topology->angles[i].atom3;

    const Vector3D &posi1 = app->positions[a1];
    const Vector3D &posi2 = app->positions[a2];
    const Vector3D &posi3 = app->positions[a3];
    const Vector3D &vel1 = app->velocities[a1];
    const Vector3D &vel2 = app->velocities[a2];
    const Vector3D &vel3 = app->velocities[a3];

    //now handle atoms 1 and 2
    Vector3D unitVec = posi2 - posi1;   // now it is only a vector
    unitVec.normalize();
    // now the vector is unit vector, the dist is the length of the original
    // vector. We assume the weight, w1, and w2 to be 1 always, thus the 
    // following still holds true
    // w1 = w2^2;

    Real coeff = -myGamma * ((vel2 - vel1).dot(unitVec));
    (*myDissipativeForces)[a2] += unitVec * coeff;
    (*myDissipativeForces)[a1] -= unitVec * coeff;
    Real randNum = randomGaussian(sdv, mySeed);
    coeff = mySigma * randNum;
    (*myRandomForces)[a2] += unitVec * coeff;
    (*myRandomForces)[a1] -= unitVec * coeff;

    //now handle atoms 2 and 3
    unitVec = posi3 - posi2;   // now it is only a vector
    unitVec.normalize();

    coeff = -myGamma * ((vel3 - vel2).dot(unitVec));
    (*myDissipativeForces)[a3] += unitVec * coeff;
    (*myDissipativeForces)[a2] -= unitVec * coeff;
    randNum = randomGaussian(sdv, mySeed);
    coeff = mySigma * randNum;
    (*myRandomForces)[a3] += unitVec * coeff;
    (*myRandomForces)[a2] -= unitVec * coeff;

    //now handle atoms 1 and 3
    unitVec = posi3 - posi1;   // now it is only a vector
    unitVec.normalize();

    coeff = -myGamma * ((vel3 - vel1).dot(unitVec));
    (*myDissipativeForces)[a3] += unitVec * coeff;
    (*myDissipativeForces)[a1] -= unitVec * coeff;
    randNum = randomGaussian(sdv, mySeed);
    coeff = mySigma * randNum;
    (*myRandomForces)[a3] += unitVec * coeff;
    (*myRandomForces)[a1] -= unitVec * coeff;
  }
}

void DMDLeapfrogIntegrator::run(int numTimesteps) {
  for (int i = 0; i < numTimesteps; i++) {
    preStepModify();
    doHalfKick();
    doDriftOrNextIntegrator();
    calculateForces();   // only calculate the conservative forces.
    calculateDissipativeAndRandomForces();
    doHalfKickVhat();
    //iteratively calculate the dissipative force and velocities.
    doHalfKickIterate();
    postStepModify();
  }
}

STSIntegrator *DMDLeapfrogIntegrator::doMake(const vector<Value> &values,
                                             ForceGroup *fg) const {
  return new DMDLeapfrogIntegrator(values[0], values[1], values[2], values[3],
                                   values[4], fg);
}

void DMDLeapfrogIntegrator::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("iterations",
               Value(myNumIter, ConstraintValueType::NotNegative()), 2));
  parameters.push_back
    (Parameter("gamma",
               Value(myGamma * 1000, ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("temperature",
               Value(myTemperature, ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("seed", Value(mySeed, ConstraintValueType::Positive()), 1234));
}

