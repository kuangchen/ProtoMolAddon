#include <protomol/output/OutputXYZTrajectoryForce.h>
#include <protomol/module/MainModule.h>
#include <protomol/config/Configuration.h>
#include <protomol/integrator/Integrator.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/io/XYZTrajectoryWriter.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


const string OutputXYZTrajectoryForce::keyword("XYZForceFile");


OutputXYZTrajectoryForce::OutputXYZTrajectoryForce() : xYZ() {}


OutputXYZTrajectoryForce::OutputXYZTrajectoryForce(const string &filename,
                                                   int freq) :
  Output(freq), xYZ(new XYZTrajectoryWriter(filename)) {}


OutputXYZTrajectoryForce::~OutputXYZTrajectoryForce() {
  if (xYZ) delete xYZ;
}


void OutputXYZTrajectoryForce::doInitialize() {
  if (xYZ == NULL || !xYZ->open())
    THROWS("Can not open '" << (xYZ ? xYZ->getFilename() : "")
           << "' for " << getId() << ".");
}


void OutputXYZTrajectoryForce::doRun(int) {
  if (!xYZ->write(*(app->integrator->getForces()), app->topology->atoms,
                    app->topology->atomTypes))
    THROWS("Could not write " << getId() << " '" << xYZ->getFilename()
           << "'.");
}


void OutputXYZTrajectoryForce::doFinalize(int) {
  xYZ->close();
}


Output *OutputXYZTrajectoryForce::doMake(const vector<Value> &values) const {
  return new OutputXYZTrajectoryForce(values[0], values[1]);
}


void OutputXYZTrajectoryForce::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(xYZ ? xYZ->getFilename() : "",
                              ConstraintValueType::NotEmpty())));
  Output::getParameters(parameter);
}
