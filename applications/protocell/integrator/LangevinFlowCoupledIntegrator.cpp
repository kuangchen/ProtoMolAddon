#include "LangevinFlowCoupledIntegrator.h"
#include <protomol/base/Report.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/force/ForceGroup.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/topology/TopologyUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/PMConstants.h>

using namespace std; 
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____ LangevinFlowCoupledIntegrator

const string LangevinFlowCoupledIntegrator::keyword("LangevinFlowCoupled");

LangevinFlowCoupledIntegrator::LangevinFlowCoupledIntegrator() :
  STSIntegrator(), myLangevinTemperature(-1.0), myGamma(-1.0),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(-1),
  averageVelocityX(0.0), averageVelocityY(0.0), averageVelocityZ(0.0),
  slopeVelocityX(0.0), slopeVelocityY(0.0), slopeVelocityZ(0.0)
   {}

LangevinFlowCoupledIntegrator::
LangevinFlowCoupledIntegrator(Real timestep, Real LangevinTemperature, Real gamma,
                          int seed, Real avVX, Real avVY, Real avVZ,
                          Real slVX, Real slVY, Real slVZ,
                          ForceGroup *overloadedForces) :
  STSIntegrator(timestep, overloadedForces),
  myLangevinTemperature(LangevinTemperature),
  myGamma(gamma),// / (1000 * Constant::INV_TIMEFACTOR)),
  // gamma is in Kcal/ps, myGamma is in Kcal/(fs*INV_TIMEFACTOR)
  mySeed(seed),
  averageVelocityX(avVX), averageVelocityY(avVY), averageVelocityZ(avVZ),
  slopeVelocityX(slVX), slopeVelocityY(slVY), slopeVelocityZ(slVZ)
  {}

void LangevinFlowCoupledIntegrator::initialize(ProtoMolApp *app) {

  STSIntegrator::initialize(app);
  initializeForces();

  // Get the cell centers

  const GenericTopology *topo = app->topology;
  const unsigned int csize = topo->atoms.size();

  //loop to find cells
  for (unsigned int i=0; i < csize; i++) {

    //by convention CA, GA etc.
    if((topo->atoms[i].name.c_str())[1] == 'A'){

      //report << hint << "Center Atom name " << i << endr;
      //save center atom as index by residue
      cellCenters[topo->atoms[i].residue_seq] = i;

    }
    
  }

  //get size=number of cells
  numCells = cellCenters.size();

  //loop to test for cell centers cells
  for (unsigned int i=0; i < csize; i++)
      cellCenters[topo->atoms[i].residue_seq];

  //if more cels than centers size will be greater as map will be extended
  if( cellCenters.size() != numCells || cellCenters.size() == 0 ){
      report << error << "Not all cell centers defined!"<<endr;
  }

}

void LangevinFlowCoupledIntegrator::doDrift() {
  const Real h = getTimestep();// * Constant::INV_TIMEFACTOR;
  app->positions.intoWeightedAdd(h, app->velocities);
  buildMolecularCenterOfMass(&app->positions, app->topology);
  buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinFlowCoupledIntegrator::doHalfKick() {
    const unsigned int count = app->positions.size();
    const Real dt = getTimestep();// * Constant::INV_TIMEFACTOR; // in fs

    //variables for diagnostics
    Real average_velocity = 0.0;
    Real ke = 0.0;

    for (unsigned int i = 0; i < count; i++ ) {

        //user defined fixed velocity
        Vector3D fluidVelocity(
                               slopeVelocityX * app->positions[i][0] + averageVelocityX, 
                               slopeVelocityY * app->positions[i][1] + averageVelocityY, 
                               slopeVelocityZ * app->positions[i][2] + averageVelocityZ);

        //find damping factor relative to the
        //vector connecting thr SCE to the cell center
        const unsigned int myCellCenter = cellCenters[app->topology->atoms[i].residue_seq];
        const Real normFluidV = fluidVelocity.norm();
        Real factor = 1.0;

        //dont use velocity if center
        if(myCellCenter != i && normFluidV != 0.0){
          Vector3D sc =  app->positions[myCellCenter] - app->positions[i];
                        //app->topology->minimalDifference( app->positions[i], app->positions[myCellCenter]);
          factor = sc.dot(fluidVelocity) / ( normFluidV * sc.norm() );
          //if(factor > 1.0) report << hint << "factor too big " << factor << endr;
          if(factor < 0.0){
            factor = 0.1;
          }
        }

        //if cell center
        if( myCellCenter == i ) factor = 0.1;

        //factor must be +ve here, so scale gamma
        Real aGamma = myGamma;// * factor;  //####Removed factor

        //Langevin leapfrog from here
        const Real fdt = ( 1.0 - exp( -0.5 * aGamma * dt ) ) / aGamma;
        const Real vdt = exp(-0.5*aGamma*dt);
        const Real ndt = sqrt( ( 1.0 - exp( -aGamma * dt ) ) / (2.0 * aGamma) );
        const Real forceConstant = 2 * Constant::SI::BOLTZMANN * 1.0e15 //for SI units, BOLTZMANN //
                                        * myLangevinTemperature 
                                            * aGamma; 

        //  Generate gaussian random numbers for each spatial direction
        Vector3D gaussRandCoord1(randomGaussianNumber(mySeed),
                                 randomGaussianNumber(mySeed),
                                 randomGaussianNumber(mySeed));
        Real mass = app->topology->atoms[i].scaledMass;
        Real sqrtFCoverM = sqrt(forceConstant / mass);

        //remove velocity
        app->velocities[i] -= fluidVelocity;

        // semi-update velocities
        app->velocities[i] = app->velocities[i]*vdt
                                +(*myForces)[i] * fdt / mass
                                  //+ projectedVelocity * fdt
                                    +gaussRandCoord1*sqrtFCoverM*ndt;
                                    //+gaussRandCoord1*(sqrt(sqrtFCoverM + sqrtFCoverMf))*ndt;

        //find "real" temperature
        for(int k=0; k<3; k++){
          ke += 0.5 * app->velocities[i].c[k] * app->velocities[i].c[k] * mass;
        }

        //replace velocity
        app->velocities[i] += fluidVelocity;

        //find average for diagnostics
        average_velocity += app->velocities[i].c[0];

    }

    //diagnostic output
    //report << hint << "Average velocity  " << average_velocity / (Real)count
    //        << " Set velocity " << averageVelocityX << " Temp " << 2.0 * ke / Constant::BOLTZMANN / count / 3.0 << endr;

    buildMolecularMomentum(&app->velocities, app->topology);
}

void LangevinFlowCoupledIntegrator::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (Parameter("temperature", Value(myLangevinTemperature,
                                    ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("gamma", Value(myGamma /** (1000 * Constant::INV_TIMEFACTOR)*/,
                              ConstraintValueType::NotNegative())));
  parameters.push_back
    (Parameter("seed", Value(mySeed, ConstraintValueType::NotNegative()),
               1234));
  parameters.push_back
    (Parameter("averageflowvelocityx", Value(averageVelocityX, ConstraintValueType::NoConstraints()),
             0.0, Text("Average flow velocity x")));
  parameters.push_back
    (Parameter("averageflowvelocityy", Value(averageVelocityY, ConstraintValueType::NoConstraints()),
             0.0, Text("Average flow velocity y")));
  parameters.push_back
    (Parameter("averageflowvelocityz", Value(averageVelocityZ, ConstraintValueType::NoConstraints()),
             0.0, Text("Average flow velocity z")));
  parameters.push_back
    (Parameter("slopeflowvelocityx", Value(slopeVelocityX, ConstraintValueType::NoConstraints()),
             0.0, Text("Slope of flow velocity x")));
  parameters.push_back
    (Parameter("slopeflowvelocityy", Value(slopeVelocityY, ConstraintValueType::NoConstraints()),
             0.0, Text("Slope of flow velocity y")));
  parameters.push_back
    (Parameter("slopeflowvelocityz", Value(slopeVelocityZ, ConstraintValueType::NoConstraints()),
             0.0, Text("Slope of flow velocity z")));
}

STSIntegrator *LangevinFlowCoupledIntegrator::doMake(const vector<Value> &values,
                                                 ForceGroup *fg) const {
  return new LangevinFlowCoupledIntegrator(values[0], values[1], values[2],
                                            values[3], values[4], values[5], values[6],
                                            values[7], values[8], values[9],
                                            fg);
}
