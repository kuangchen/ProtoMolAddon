#include <protomol/output/OutputXYZTrajectoryPos.h>
#include <protomol/module/MainModule.h>
#include <protomol/output/OutputCache.h>
#include <protomol/config/Configuration.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/base/Exception.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/io/XYZTrajectoryWriter.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


const string OutputXYZTrajectoryPos::keyword("XYZPosFile");

OutputXYZTrajectoryPos::OutputXYZTrajectoryPos() :
  xYZ(0), minimalImage(false) {}


OutputXYZTrajectoryPos::OutputXYZTrajectoryPos(const string &filename, int freq,
                                               bool minimal) :
  Output(freq), xYZ(new XYZTrajectoryWriter(filename)),
  minimalImage(minimal) {}


OutputXYZTrajectoryPos::~OutputXYZTrajectoryPos() {
  if (xYZ) delete xYZ;
}


void OutputXYZTrajectoryPos::doInitialize() {
  if (!xYZ || !xYZ->open())
    THROWS("Can not open '" << (xYZ ? xYZ->getFilename() : "")
           << "' for " << getId() << ".");
}


void OutputXYZTrajectoryPos::doRun(int) {
  const Vector3DBlock *pos =
    (minimalImage ? app->outputCache.getMinimalPositions() : &app->positions);

  if (!xYZ->write(*pos, app->topology->atoms, app->topology->atomTypes))
    THROWS("Could not write " << getId() << " '" << xYZ->getFilename()
           << "'.");
}


void OutputXYZTrajectoryPos::doFinalize(int) {
  xYZ->close();
}


Output *OutputXYZTrajectoryPos::doMake(const vector<Value> &values) const {
  return new OutputXYZTrajectoryPos(values[0], values[1], values[2]);
}


void OutputXYZTrajectoryPos::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(xYZ != NULL ? xYZ->getFilename() : "",
                              ConstraintValueType::NotEmpty())));

  Output::getParameters(parameter);

  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(minimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not")));
}


bool OutputXYZTrajectoryPos::
adjustWithDefaultParameters(vector<Value> &values,
                            const Configuration *config) const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
    values[2] = (*config)[InputMinimalImage::keyword];

  return checkParameters(values);
}
