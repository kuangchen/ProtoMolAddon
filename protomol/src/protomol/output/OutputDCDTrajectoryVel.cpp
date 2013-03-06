#include <protomol/output/OutputDCDTrajectoryVel.h>
#include <protomol/config/Configuration.h>
#include <protomol/output/OutputCache.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/io/DCDTrajectoryWriter.h>
#include <protomol/base/Exception.h>
#include <protomol/ProtoMolApp.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


const string OutputDCDTrajectoryVel::keyword("DCDVELFile");


OutputDCDTrajectoryVel::OutputDCDTrajectoryVel() :
  dCD(0), minimalImage(false) {}


OutputDCDTrajectoryVel::OutputDCDTrajectoryVel(const string &filename,
                                               int freq, bool minimal) :
  Output(freq), dCD(new DCDTrajectoryWriter(filename)),
  minimalImage(minimal) {}


OutputDCDTrajectoryVel::~OutputDCDTrajectoryVel() {
  if (dCD) delete dCD;
}


void OutputDCDTrajectoryVel::doInitialize() {
  if (!dCD || !dCD->open())
    THROWS("Can not open '" << (dCD ? dCD->getFilename() : "")
           << "' for " << getId() << ".");
}

void OutputDCDTrajectoryVel::doRun(int) {
  if (!dCD->write(app->velocities))
    THROWS("Could not write " << getId() << " '" << dCD->getFilename()
           << "'.");
}


void OutputDCDTrajectoryVel::doFinalize(int) {
  dCD->close();
}


Output *OutputDCDTrajectoryVel::doMake(const vector<Value> &values) const {
  return new OutputDCDTrajectoryVel(values[0], values[1], values[2]);
}


void OutputDCDTrajectoryVel::getParameters(vector<Parameter> &parameter)
const {
  parameter.push_back
    (Parameter(getId(), Value(dCD ? dCD->getFilename() : "",
                              ConstraintValueType::NotEmpty())));
  Output::getParameters(parameter);
  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(minimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not (NA for Vel)")));
}

bool OutputDCDTrajectoryVel::adjustWithDefaultParameters(
  vector<Value> &values, const Configuration *config) const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
    values[2] = (*config)[InputMinimalImage::keyword];

  return checkParameters(values);
}
