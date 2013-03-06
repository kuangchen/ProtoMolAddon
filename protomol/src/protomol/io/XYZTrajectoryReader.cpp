#include <protomol/io/XYZTrajectoryReader.h>

#include <protomol/type/String.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol;

//____XYZTrajectoryReader

XYZTrajectoryReader::XYZTrajectoryReader() : xyz(0), first(true) {}

XYZTrajectoryReader::XYZTrajectoryReader(const string &filename) :
  XYZReader(filename), xyz(0), first(true) {}

XYZTrajectoryReader::~XYZTrajectoryReader() {
  if (xyz) delete xyz;
}

bool XYZTrajectoryReader::tryFormat() {
  Vector3DBlock xyz;
  return read(xyz);
}

bool XYZTrajectoryReader::read() {
  if (!xyz) xyz = new Vector3DBlock();
  return read(*xyz);
}

bool XYZTrajectoryReader::read(Vector3DBlock &xyz) {
  try {
    doRead(xyz);
    return true;

  } catch (const Exception &e) {}

  return false;
}

void XYZTrajectoryReader::doRead(Vector3DBlock &xyz) {
  if (!is_open()) {
    if (!open()) THROW("open failed");

  } else file.seekg(0);

  vector<string> tokens;

  // Number of frames
  if (first)
    if (getLineTokens(tokens) != 1) THROW("Invalid frame count");

  try {
    // Read frames
    //for (unsigned int i = 0; i < n && !file.fail(); ++i)
    vector<string> names;
    XYZReader::doRead(xyz, names);
    
    if (file.fail()) THROW("Reading data failed");

  } catch (const Exception &e) {
    xyz.clear();
    throw e;
  }

  first = false;
}

Vector3DBlock *XYZTrajectoryReader::orphanXYZ() {
  Vector3DBlock *tmp = xyz;
  xyz = 0;
  return tmp;
}

namespace ProtoMol {
  XYZTrajectoryReader &operator>>(XYZTrajectoryReader &reader,
                                  Vector3DBlock &xyz) {
    reader.doRead(xyz);
    return reader;
  }
}
