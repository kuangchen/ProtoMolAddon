#include <protomol/output/OutputFinalXYZVel.h>
#include <protomol/config/Configuration.h>
#include <protomol/output/OutputCache.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/io/XYZWriter.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


const string OutputFinalXYZVel::keyword("finXYZVelFile");


OutputFinalXYZVel::OutputFinalXYZVel() : Output(-1) {}


OutputFinalXYZVel::OutputFinalXYZVel(const string &filename) :
  Output(-1), filename(filename) {}


void OutputFinalXYZVel::doFinalize(int step) {
  XYZWriter writer;
  if (!writer.open(filename))
    THROWS("Can't open " << getId() << " '" << filename << "'.");

  writer.setComment("Time : " + toString(app->outputCache.getTime()) +
                    ", step : " + toString(step) + ".");

  if (!writer.write(*&app->velocities, app->topology->atoms,
                    app->topology->atomTypes))
    THROWS("Could not write " << getId() << " '" << filename << "'.");
}


Output *OutputFinalXYZVel::doMake(const vector<Value> &values) const {
  return new OutputFinalXYZVel(values[0]);
}


void OutputFinalXYZVel::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(filename, ConstraintValueType::NotEmpty())));
}

