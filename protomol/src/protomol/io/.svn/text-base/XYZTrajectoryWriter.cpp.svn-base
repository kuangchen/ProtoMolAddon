#include <protomol/io/XYZTrajectoryWriter.h>

#include <iomanip>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____XYZTrajectoryWriter

XYZTrajectoryWriter::XYZTrajectoryWriter() :
  Writer(), myCoords(NULL), myNames(NULL), myAtoms(NULL), myAtomTypes(NULL),
  myFirst(true) {}

XYZTrajectoryWriter::XYZTrajectoryWriter(const string &filename) :
  Writer(filename), myCoords(NULL), myNames(NULL), myAtoms(NULL),
  myAtomTypes(NULL), myFirst(true) {}

XYZTrajectoryWriter::~XYZTrajectoryWriter() {}

bool XYZTrajectoryWriter::openWith(const XYZ &xyz) {
  setNames(xyz.names);
  myFirst = false;
  return open();
}

bool XYZTrajectoryWriter::openWith(const vector<string> &names) {
  setNames(names);
  myFirst = false;
  return open();
}

bool XYZTrajectoryWriter::openWith(const vector<Atom> &atoms,
                                   const vector<AtomType> &atomTypes) {
  setNames(atoms, atomTypes);
  myFirst = false;
  return open();
}

bool XYZTrajectoryWriter::openWith(const string &filename, const XYZ &xyz) {
  setNames(xyz.names);
  myFirst = false;
  return open(filename);
}

bool XYZTrajectoryWriter::openWith(const string &filename,
                                   const vector<string> &names) {
  setNames(names);
  myFirst = false;
  return open(filename);
}

bool XYZTrajectoryWriter::openWith(const string &filename,
                                   const vector<Atom> &atoms,
                                   const vector<AtomType> &atomTypes) {
  setNames(atoms, atomTypes);
  myFirst = false;
  return open(filename);
}

bool XYZTrajectoryWriter::reopen() {
  if (myFirst) {
    myFirst = false;
    if (!open()) return false;
  }

  if (is_open()) close();

  // Try to read the number of frames
  file.clear();
  open(filename.c_str(), ios::in);
  string line, str;
  line = getline();
  close();
  stringstream ss(line);
  ss >> str;
  if (toInt(str) > 0) {
    // Ok, we have already written frames
    str = (toString(toInt(str) + 1) + "                       ").substr(0, 19);
    str += "\n";
    file.clear();
    open(filename.c_str(), ios::in | ios::out);
    file.seekp(0, ios::beg);
    file.write(str.c_str(), str.size());

  } else {
    // First time ...
    file.clear();
    open(filename.c_str(), ios::out | ios::trunc);
    str = (toString(1) + "                       ").substr(0, 19);
    str += "\n";
    file << str;
  }
  close();
  file.clear();
  open(filename.c_str(), ios::out | ios::app);

  return !file.fail();
}

bool XYZTrajectoryWriter::write(const XYZ &xyz) {
  setCoords(xyz.coords);
  setNames(xyz.names);
  return write();
}

bool XYZTrajectoryWriter::write(const Vector3DBlock &coords) {
  setCoords(coords);
  return write();
}

bool XYZTrajectoryWriter::write(const Vector3DBlock &coords,
                                const vector<string> &names) {
  setNames(names);
  setCoords(coords);
  return write();
}

bool XYZTrajectoryWriter::write(const Vector3DBlock &coords,
                                const vector<Atom> &atoms,
                                const vector<AtomType> &atomTypes) {
  setNames(atoms, atomTypes);
  setCoords(coords);
  return write();
}

bool XYZTrajectoryWriter::write() {
  if (!reopen()) return false;

  if (myCoords == NULL)
    report << error << "[XYZTrajectoryWriter::write]"
           << " No coordinates specified." << endr;

  if (!(myNames != NULL || (myAtoms != NULL && myAtomTypes != NULL)))
    report << error << "[XYZTrajectoryWriter::write]"
           << " No atom names specified." << endr;

  const unsigned int count = myCoords->size();
  if ((myNames != NULL && myNames->size() != count) ||
      (myAtoms != NULL && myAtomTypes != NULL && myAtoms->size() != count))
    report << error << "[XYZTrajectoryWriter::write]"
           << " Coorindate and atom name size are not equal." << endr;


  // First, write the number of atoms
  file << count << endl;

  // Write atoms
  file << setprecision(15);   // This should be some FLT_DIG or DBL_DIG ...
  for (unsigned int i = 0; i < count; ++i) {
    file << (myNames ? (*myNames)[i] : (*myAtomTypes)[(*myAtoms)[i].type].name)
         << "\t";
    file.width(24);
    file << (*myCoords)[i].c[0];
    file.width(24);
    file << (*myCoords)[i].c[1];
    file.width(24);
    file << (*myCoords)[i].c[2];
    file << endl;
  }

  close();
  return !file.fail();
}

void XYZTrajectoryWriter::setNames(const XYZ &xyz) {
  myNames = &(xyz.names);
  myAtoms = NULL;
  myAtomTypes = NULL;
}

void XYZTrajectoryWriter::setNames(const vector<string> &names) {
  myNames = &names;
  myAtoms = NULL;
  myAtomTypes = NULL;
}

void XYZTrajectoryWriter::setNames(const vector<Atom> &atoms,
                                   const vector<AtomType> &atomTypes) {
  myNames = NULL;
  myAtoms = &atoms;
  myAtomTypes = &atomTypes;
}

void XYZTrajectoryWriter::setCoords(const Vector3DBlock &coords) {
  myCoords = &coords;
}

namespace ProtoMol {
  XYZTrajectoryWriter &operator<<(XYZTrajectoryWriter &xyzWriter,
                                  const XYZ &xyz) {
    xyzWriter.write(xyz);
    return xyzWriter;
  }
  
  XYZTrajectoryWriter &operator<<(XYZTrajectoryWriter &xyzWriter,
                                  const Vector3DBlock &coords) {
    xyzWriter.write(coords);
    return xyzWriter;
  }
}
