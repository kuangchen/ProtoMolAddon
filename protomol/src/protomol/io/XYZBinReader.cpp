#include <protomol/io/XYZBinReader.h>

#include <protomol/base/Report.h>
#include <protomol/base/SystemUtilities.h>
#include <protomol/type/TypeSelection.h>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;
//____XYZBinReader

XYZBinReader::XYZBinReader() :
  Reader(ios::binary), myCoords(NULL) {}

XYZBinReader::XYZBinReader(const string &filename) :
  Reader(ios::binary, filename), myCoords(NULL) {}

XYZBinReader::~XYZBinReader() {
  if (myCoords != NULL) delete myCoords;
}

bool XYZBinReader::tryFormat() {
  if (!open()) return false;

  file.seekg(0, ios::end);
  ios::pos_type size = file.tellg();
  size -= 4;
  file.seekg(0, ios::beg);

  typedef TypeSelection::Int<4>::type int32;
  int32 n = 0;
  File::read(reinterpret_cast<char *>(&n), sizeof(int32));
  close();

  if (static_cast<int32>(size / (3 * sizeof(float))) == n ||
      static_cast<int32>(size / (3 * sizeof(double))) == n ||
      static_cast<int32>(size / (3 * sizeof(Real))) == n)
    return !file.fail();
  else {
    swapBytes(n);
    if (static_cast<int32>(size / (3 * sizeof(float))) == n ||
        static_cast<int32>(size / (3 * sizeof(double))) == n ||
        static_cast<int32>(size / (3 * sizeof(Real))) == n)
      return !file.fail();
  }

  return false;
}

bool XYZBinReader::read() {
  if (myCoords == NULL)
    myCoords = new Vector3DBlock();
  return read(*myCoords);
}

bool XYZBinReader::read(Vector3DBlock &coords) {
  if (!tryFormat()) return false;
  if (!open()) return false;

  file.seekg(0, ios::end);
  ios::pos_type size = file.tellg();
  size -= 4;
  file.seekg(0, ios::beg);
  bool swapEndian = false;

  typedef TypeSelection::Int<4>::type int32;
  int32 n = 0;
  File::read(reinterpret_cast<char *>(&n), sizeof(int32));

  if (static_cast<int32>(size / (3 * sizeof(float))) != n &&
      static_cast<int32>(size / (3 * sizeof(double))) != n &&
      static_cast<int32>(size / (3 * sizeof(Real))) != n) {
    swapBytes(n);
    swapEndian = true;
    if (static_cast<int32>(size / (3 * sizeof(float))) == n ||
        static_cast<int32>(size / (3 * sizeof(double))) == n ||
        static_cast<int32>(size / (3 * sizeof(Real))) == n)
      report << hint << "[XYZBinReader::read] Reading " <<
      (ISLITTLEENDIAN ? "big" : "little")
             << "endian input on " << (ISLITTLEENDIAN ? "little" : "big") <<
      "endian machine." << endr;
  }

  if (static_cast<int32>(size / (3 * sizeof(Real))) == n) {
    coords.resize(n);

    Real *vec = new Real[n * 3];
    File::read(reinterpret_cast<char *>(vec), n * 3 * sizeof(Real));

    for (int32 i = 0; i < n; ++i) {
      if (swapEndian) {
        swapBytes(vec[i * 3 + 0]);
        swapBytes(vec[i * 3 + 1]);
        swapBytes(vec[i * 3 + 2]);
      }
      coords[i].c[0] = static_cast<Real>(vec[i * 3 + 0]);
      coords[i].c[1] = static_cast<Real>(vec[i * 3 + 1]);
      coords[i].c[2] = static_cast<Real>(vec[i * 3 + 2]);
    }

    delete[] vec;
  } else if (static_cast<int32>(size / (3 * sizeof(float))) == n) {
    coords.resize(n);

    float *vec = new float[n * 3];
    File::read(reinterpret_cast<char *>(vec), n * 3 * sizeof(float));

    for (int32 i = 0; i < n; ++i) {
      if (swapEndian) {
        swapBytes(vec[i * 3 + 0]);
        swapBytes(vec[i * 3 + 1]);
        swapBytes(vec[i * 3 + 2]);
      }
      coords[i].c[0] = static_cast<Real>(vec[i * 3 + 0]);
      coords[i].c[1] = static_cast<Real>(vec[i * 3 + 1]);
      coords[i].c[2] = static_cast<Real>(vec[i * 3 + 2]);
    }

    delete[] vec;
    report << hint << "[XYZBinReader::read] Conversion from " <<
    sizeof(float) << " to " << sizeof(Real) << " bytes." << endr;
  } else if (static_cast<int32>(size / (3 * sizeof(double))) == n) {
    coords.resize(n);

    double *vec = new double[n * 3];
    File::read(reinterpret_cast<char *>(vec), n * 3 * sizeof(double));

    for (int32 i = 0; i < n; ++i) {
      if (swapEndian) {
        swapBytes(vec[i * 3 + 0]);
        swapBytes(vec[i * 3 + 1]);
        swapBytes(vec[i * 3 + 2]);
      }
      coords[i].c[0] = static_cast<Real>(vec[i * 3 + 0]);
      coords[i].c[1] = static_cast<Real>(vec[i * 3 + 1]);
      coords[i].c[2] = static_cast<Real>(vec[i * 3 + 2]);
    }

    delete[] vec;
    report << hint << "[XYZBinReader::read] Conversion from " <<
    sizeof(double) << " to " << sizeof(Real) << " bytes." << endr;
  } else {
    report << recoverable << "[XYZBinReader::read]"
           << " XYZBin file \'" << filename
           << "\' could find adequate float nor double type." << endr;
    close();
    return false;
  }

  cout << "Success reading XYZBin" << endl;

  close();
  return !file.fail();
}

XYZ XYZBinReader::getXYZ() const {
  XYZ res;
  if (myCoords != NULL)
    res.coords = (*myCoords);
  res.names.resize(res.coords.size(), "NONAME");
  return res;
}

Vector3DBlock *XYZBinReader::orphanCoords() {
  Vector3DBlock *tmp = myCoords;
  myCoords = NULL;
  return tmp;
}

namespace ProtoMol {
  XYZBinReader &operator>>(XYZBinReader &xyzbinReader, XYZ &xyz) {
    xyzbinReader.read(xyz.coords);
    if (xyz.coords.size() != xyz.names.size())
      xyz.names.resize(xyz.coords.size(), "NONAME");
    return xyzbinReader;
  }

  XYZBinReader &operator>>(XYZBinReader &xyzbinReader,
                           Vector3DBlock &coords) {
    xyzbinReader.read(coords);
    return xyzbinReader;
  }
}

