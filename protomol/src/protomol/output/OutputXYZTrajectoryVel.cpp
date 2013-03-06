#include <protomol/output/OutputXYZTrajectoryVel.h>
#include <protomol/config/Configuration.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/Exception.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/io/XYZTrajectoryWriter.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


const string OutputXYZTrajectoryVel::keyword("XYZVelFile");


OutputXYZTrajectoryVel::OutputXYZTrajectoryVel() : xYZ() {}


OutputXYZTrajectoryVel::OutputXYZTrajectoryVel(const string &filename,
                                               int freq) :
  Output(freq), xYZ(new XYZTrajectoryWriter(filename)) {}


OutputXYZTrajectoryVel::~OutputXYZTrajectoryVel() {
  if (xYZ) delete xYZ;
}


void OutputXYZTrajectoryVel::doInitialize() {
  if (!xYZ || !xYZ->open())
    THROWS("Can not open '" << (xYZ ? xYZ->getFilename() : "")
           << "' for " << getId() << ".");
}


void OutputXYZTrajectoryVel::doRun(int) {
  if (!xYZ->write(*&app->velocities, app->topology->atoms,
                    app->topology->atomTypes))
    THROWS("Could not write " << getId() << " '" << xYZ->getFilename()
           << "'.");
}


void OutputXYZTrajectoryVel::doFinalize(int) {
  xYZ->close();
}


Output *OutputXYZTrajectoryVel::doMake(const vector<Value> &values) const {
  return new OutputXYZTrajectoryVel(values[0], values[1]);
}


void OutputXYZTrajectoryVel::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(xYZ != NULL ? xYZ->getFilename() : "",
                              ConstraintValueType::NotEmpty())));

  Output::getParameters(parameter);
}
