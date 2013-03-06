#include <protomol/io/XYZWriter.h>

#include <iomanip>

#include <protomol/base/Report.h>
#include <protomol/base/StringUtilities.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____XYZWriter

XYZWriter::XYZWriter() :
  Writer() {}

XYZWriter::XYZWriter(const string &filename) :
  Writer(filename) {}

bool XYZWriter::write(const XYZ &xyz) {
  return write(xyz.coords, xyz.names);
}

bool XYZWriter::write(const Vector3DBlock &coords, const vector<Atom> &atoms,
                      const vector<AtomType> &atomTypes) {
  const unsigned int count = atoms.size();
  vector<string> names(count);
  for (unsigned int i = 0; i < count; ++i)
    names[i] = atomTypes[atoms[i].type].name;

  return write(coords, names);
}

bool XYZWriter::write(const Vector3DBlock &coords,
                      const vector<string> &names) {
  if (!open())
    return false;

  const unsigned int count = coords.size();
  if (names.size() != count)
    report << error << "[XYZWriter::write]"
           << " Coorindate and atom name size are not equal." << endr;

  // First, write the number of atoms
  file << count << endl;

  // Write atoms
  file.precision(15);// << setprecision(15);   // This should be some FLT_DIG or DBL_DIG ...
  for (unsigned int i = 0; i < count; ++i) {
    file << names[i] << "\t";
    file.width(24);
    file << coords[i].c[0];
    file.width(24);
    file << coords[i].c[1];
    file.width(24);
    file << coords[i].c[2];
    file << endl;
  }

  close();
  return !file.fail();
}


namespace ProtoMol {
  XYZWriter &operator<<(XYZWriter &xyzWriter, const XYZ &xyz) {
    xyzWriter.write(xyz.coords, xyz.names);
    return xyzWriter;
  }
}
