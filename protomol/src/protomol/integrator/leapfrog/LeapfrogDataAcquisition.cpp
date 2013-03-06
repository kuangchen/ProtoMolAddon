#include <protomol/integrator/leapfrog/LeapfrogDataAcquisition.h>
#include <protomol/io/DCDTrajectoryWriter.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;

const string LeapfrogDataAcquisition::keyword("LeapfrogDataAcquisition");

LeapfrogDataAcquisition::LeapfrogDataAcquisition() :
  LeapfrogIntegrator(), numNonWaters(0) {
    myWriterF = myWriterX = myWriterV = NULL;
    nonWaterF = nonWaterX = nonWaterV = NULL;
  }

LeapfrogDataAcquisition::LeapfrogDataAcquisition(Real timestep,
                     string dcdfile, 
                     ForceGroup *overloadedForces) :
  LeapfrogIntegrator(timestep, overloadedForces), myDCDFile(dcdfile), numNonWaters(0) {

    nonWaterF = nonWaterX = nonWaterV = NULL;

    int dcdpos = dcdfile.rfind(".dcd");
    string dcdx = dcdfile; dcdx.insert(dcdpos,"Pos");
    string dcdf = dcdfile; dcdf.insert(dcdpos,"Force");
    string dcdv = dcdfile; dcdv.insert(dcdpos,"Vel");

    myWriterF = new DCDTrajectoryWriter(dcdf, timestep, 0, ISLITTLEENDIAN);
    myWriterX = new DCDTrajectoryWriter(dcdx, timestep, 0, ISLITTLEENDIAN);
    myWriterV = new DCDTrajectoryWriter(dcdv, timestep, 0, ISLITTLEENDIAN);
}

LeapfrogDataAcquisition::~LeapfrogDataAcquisition() {
  zap(myWriterF);  
  zap(myWriterX);
  zap(myWriterV);
  zap(nonWaterF);  
  zap(nonWaterX);
  zap(nonWaterV);
}

void LeapfrogDataAcquisition::initialize(ProtoMolApp *app) {
  STSIntegrator::initialize(app);
  initializeForces();
  for (unsigned int i = 0; i < myForces->size(); i++)
    if (!app->topology->molecules[app->topology->atoms[i].molecule].water)
      numNonWaters++;
  nonWaterF = new Vector3DBlock(numNonWaters);
  nonWaterX = new Vector3DBlock(numNonWaters);
  nonWaterV = new Vector3DBlock(numNonWaters);
}

void LeapfrogDataAcquisition::run(int numTimesteps) {

  if (numTimesteps < 1)
    return;
  preStepModify();
  doHalfKickdoDrift();
  calculateForces();
  writeDCD();
  for (int i = 1; i < numTimesteps; i++) {
    doKickdoDrift();
    calculateForces();
    writeDCD();
  }

  doHalfKick();
  postStepModify();

}

void LeapfrogDataAcquisition::writeDCD()
{
  // Create a structure to hold force
  // values for non-water molecules.
  int pos = 0;
  for (unsigned int i = 0; (i < app->topology->atoms.size()) && pos < numNonWaters; i++) 
    if (!app->topology->molecules[app->topology->atoms[i].molecule].water) {
      (*nonWaterF)[pos] = (*myForces)[i];
      (*nonWaterX)[pos] = app->positions[i];
      (*nonWaterV)[pos] = app->velocities[i];
      pos++;
    }   
  // Write the DCD.
  myWriterF->write(*nonWaterF);
  myWriterX->write(*nonWaterX);
  myWriterV->write(*nonWaterV);
}

void LeapfrogDataAcquisition::getParameters(vector<Parameter> &parameters)
const {
  STSIntegrator::getParameters(parameters);
  parameters.push_back
    (
    Parameter("dcdfile", Value(myDCDFile), ConstraintValueType::NotEmpty())
    );
}

STSIntegrator *LeapfrogDataAcquisition::doMake(const vector<Value> &values,
                       ForceGroup *fg) const {
  return new LeapfrogDataAcquisition(values[0], values[1], fg);
}
