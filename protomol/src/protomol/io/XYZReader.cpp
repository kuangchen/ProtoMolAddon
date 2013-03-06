#include <protomol/io/XYZReader.h>

#include <protomol/type/String.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol;
//____XYZReader

XYZReader::XYZReader() :
  Reader(), myCoords(NULL), myNames(NULL) {}

XYZReader::XYZReader(const string &filename) :
  Reader(filename), myCoords(NULL), myNames(NULL) {}

XYZReader::~XYZReader() {
  if (myCoords) delete myCoords;
  if (myNames) delete myNames;
}

bool XYZReader::tryFormat() {
  Vector3DBlock coords;
  vector<string> names;

  bool retVal = read(coords, names);

  close();

  return retVal;
}

bool XYZReader::read() {
  if (!myCoords) myCoords = new Vector3DBlock();
  if (!myNames) myNames = new vector<string>();
  return read(*myCoords, *myNames);
}

bool XYZReader::read(XYZ &xyz) {
  return read(xyz.coords, xyz.names);
}

bool XYZReader::read(Vector3DBlock &coords, vector<string> &names) {
  try {
    doRead(coords, names);
    return true;

  } catch (const Exception &e) {}

  return false;
}

void XYZReader::doRead(Vector3DBlock &coords, vector<string> &names) {
  if (!is_open()) if (!open()) THROW("Open failed");

  vector<string> tokens;

  // Number of atoms
  if (getLineTokens(tokens) != 1) THROW("Invalid atom count");
  unsigned int n = String::parseUInteger(tokens[0]);

  // Check for comment
  if (file.peek() == '!') comment = getline();

  coords.resize(n);
  names.resize(n);

  // Read atoms
  for (unsigned int i = 0; i < n && !file.fail(); i++) {
    unsigned int count = getLineTokens(tokens);
    if (count != 4) THROWS("Invalid XYZ line " << i << " tokens " << count);
    names[i] = tokens[0];
    coords[i].c[0] = String::parseDouble(tokens[1]);
    coords[i].c[1] = String::parseDouble(tokens[2]);
    coords[i].c[2] = String::parseDouble(tokens[3]);
  }

  if (file.fail()) THROW("Data read failed");
}

XYZ XYZReader::getXYZ() const {
  XYZ res;
  if (myCoords) res.coords = *myCoords;
  if (myNames) res.names = *myNames;
  return res;
}

Vector3DBlock *XYZReader::orphanCoords() {
  Vector3DBlock *tmp = myCoords;
  myCoords = NULL;
  return tmp;
}

vector<string> *XYZReader::orphanNames() {
  vector<string> *tmp = myNames;
  myNames = NULL;
  return tmp;
}

namespace ProtoMol {
  XYZReader &operator>>(XYZReader &reader, XYZ &xyz) {
    reader.doRead(xyz.coords, xyz.names);
    return reader;
  }
  
  XYZReader &operator>>(XYZReader &reader, Vector3DBlock &coords) {
    if (!reader.myNames) reader.myNames = new vector<string>();
    reader.read(coords, *reader.myNames);
    return reader;
  }
}
