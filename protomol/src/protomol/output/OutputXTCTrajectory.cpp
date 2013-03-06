#include "OutputXTCTrajectory.h"

#include "OutputCache.h"

#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>

#ifdef HAVE_GROMACS
extern "C" {
#include <gromacs/xtcio.h>
}
#endif

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;


const string OutputXTCTrajectory::keyword("XTCFile");


OutputXTCTrajectory::OutputXTCTrajectory() :
  fxtc(0), minimalImage(false), frameOffset(0) {}


OutputXTCTrajectory::OutputXTCTrajectory(const string &filename, int freq,
                                         bool minimal, int frameoffs) :
  Output(freq), fxtc(0), minimalImage(minimal), frameOffset(frameoffs),
  filename(filename) {

  report << plain
         << "XTC FrameOffset parameter set to " << frameOffset << "." << endr;
}


void OutputXTCTrajectory::doInitialize() {
#ifdef HAVE_GROMACS
  //  Get first frame (must exist or error)
  const int firstframe = toInt(app->config["firststep"]);

  report << debug(2) << "Firstframe " << firstframe << "." << endr;

  //  Open fil
  //  If frameOffset is zero default to overwrite data
  //  NOTE: now include "firstframe" data.
  fxtc = open_xtc(filename.c_str(), frameOffset ? "a" : "w");

  //  Test opened
  if (!fxtc) THROWS("Can not open '" << filename << "' for " << getId() << ".");

#else
  THROW("GROMACS XTC format output not available.");
#endif
}


void OutputXTCTrajectory::doRun(int step) {
#ifdef HAVE_GROMACS
  const Vector3DBlock *pos =
    (minimalImage ? app->outputCache.getMinimalPositions() : &app->positions);

  //  Bounding Box
  //  The computational box which is stored as a set of three basis vectors,
  //  to allow for triclinic PBC. For a rectangular box the box edges ar
  //  stored on the diagonal of the matrix.
  Vector3D a, b;
  app->topology->getBoundingbox(*pos, a, b);
  matrix box = {{a.c[0], a.c[1], a.c[2]}, {b.c[0], b.c[1], b.c[2]},
                {b.c[0] - a.c[0], b.c[1] - a.c[1], b.c[2] - a.c[2]}};

  //  Gromacs XYZ data struct
  const unsigned possize = pos->size(); //  Size
  SmartPointer<rvec>::Array x = new rvec[possize];

  //  Copy & convert dat
  for (unsigned i = 0; i < possize; i++)
    for (unsigned j = 0; j < 3; j++)
      x[i][j] = (*pos)[i][j] * Constant::ANGSTROM_NM;

  //  Defines precision, 1000 is the GROMACS default.
  //  Can be read from TPR file.
  real prec = 1000;
  int natoms = possize; //  Number of atom
	real time = app->outputCache.getTime() * Constant::FS_PS; //  Real time

  //  Write to fil
  if (!write_xtc((t_fileio *)fxtc, natoms, app->currentStep, time,
                 box, x.get(), prec))
    THROWS("Could not write " <<  getId() << " '" << filename << "'.");
#endif
}


void OutputXTCTrajectory::doFinalize(int) {
#ifdef HAVE_GROMACS
  close_xtc((t_fileio *)fxtc);
#endif
}


Output *OutputXTCTrajectory::doMake(const vector<Value> &values) const {
  return new OutputXTCTrajectory(values[0], values[1], values[2], values[3]);
}


void OutputXTCTrajectory::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(filename, ConstraintValueType::NotEmpty())));
  Output::getParameters(parameter);
  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(minimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not")));
  parameter.push_back
    (Parameter(keyword + "FrameOffset",
               Value(frameOffset, ConstraintValueType::NotNegative()), 0));
}


bool OutputXTCTrajectory::
adjustWithDefaultParameters(vector<Value> &values,
                            const Configuration *config) const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
    values[2] = (*config)[InputMinimalImage::keyword];

  return checkParameters(values);
}

